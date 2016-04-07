#!/usr/bin/perl
# Input a regions BED file and a BED file generated from the sequence BAM file processed through bamToBed,
# and output tables indicating the median coverage for each amplicon, each strand covered by the amplicon,
# the identities of amplicons that did not meet the minumum threshold, and a summary table for the sample.
#
# 3/4/2014 - D Sims
#############################################################################################################
use warnings;
use strict;
use autodie;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use Data::Dump;
use File::Basename;
use List::Util qw(sum);
use Cwd;

my $scriptname = basename($0);
my $version = "v1.9.7_040716";
my $description = <<"EOT";
From an Regions BED file, and a BED file generated from the sequence BAM file processed through bamToBed,
generate strand coverage information for an amplicon panel.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <regions_bed> (<BAM File> | -b <BAM BED File>)
    -b, --bambed      Start with a BAM BED file generated through BED tools
    -s, --sample      Sample name
    -i, --ion         Using an Ion Torrent BED file that has been processed.
    -r, --reads       Total number of reads in the BAM file.
    -t, --threshold   Minimum number of reads threshold (DEFAULT: 450).
    -o, --outdir      Send output to custom directory.  Default is CWD.
    -v, --version     Version information
    -h, --help        Print this help information
EOT

my $help;
my $ver_info;
my $outdir = getcwd;
my $sample_name;
my $num_reads;
my $threshold = 450;
my $ion;
my $bambed;

GetOptions( "bambed|b=s"    => \$bambed,
            "outdir|o=s"    => \$outdir,
            "sample|s=s"    => \$sample_name,
            "ion|i"         => \$ion,
            "reads|r=i"     => \$num_reads,
            "threshold|t=i" => \$threshold,
            "version|v"     => \$ver_info,
            "help|h"        => \$help )
        or print $usage;

sub help {
	printf "%s - %s\n\n%s\n\n%s\n", $scriptname, $version, $description, $usage;
	exit;
}

sub version {
	printf "%s - %s\n", $scriptname, $version;
	exit;
}

help if $help;
version if $ver_info;

my ($bamfile, $regionsbed);
if ($bambed) {
    $regionsbed = shift;
} else {
    ($regionsbed, $bamfile) = @ARGV;
}

if (! $regionsbed && ($bamfile || $bambed)) {
    die "ERROR: Insufficient arguments!\n";
    print $usage;
    exit(1);
}

if ( ! -d $outdir ) {
    print "Output directory '$outdir' not found. Creating output directory: $outdir\n";
    mkdir($outdir) || die "ERROR: Can not create $outdir: $!\n";
}

#########------------------------------ END ARG Parsing ---------------------------------#########

# If we're using an Ion Torrent processed BED file from their API, we need to process it a bit to get the Gene ID
if ( $ion ) {
    print "Converting Ion BED file to standard...\n";
    $regionsbed = proc_bed( \$regionsbed );
}

if ($bambed) {
    print "Using pre-loaded BED file: $bambed\n";
} else {
    $bambed = gen_bam_bed($bamfile);
}

# XXX
my (%coverage_data, %base_coverage_data);
get_coverage_data($bambed,$regionsbed,\%coverage_data);
get_base_coverage_data( $bambed, $regionsbed, \%base_coverage_data);

# have to count here or I get 2x counts for overlapping bases...then again, is that a bad thing?
my $total_base_reads;
$total_base_reads += $base_coverage_data{$_} for keys %base_coverage_data;

my %coverage_stats;
my @all_coverage;
get_coverage_stats(\%coverage_stats, \@all_coverage);

# XXX Generate the amplicon coverage tables
my $low_total = gen_amplicon_coverage_tables($outdir);

my ($total_bases, $total_nz_bases, $base_reads, $mean_base_coverage, $uniformity) = get_metrics(\%base_coverage_data, $total_base_reads);
# Create some output filehandles for the data
my $summary_file = "$outdir/stat_table.txt";
open( my $summary_fh, ">", "$outdir/stat_table.txt") || die "Can't open the 'stat_table.txt' file for writing: $!";

{
    # Get quartile coverage data and create output file
    select $summary_fh;
    my ( $quart1, $quart2, $quart3 ) = quartile_coverage( \@all_coverage );
    print "Sample name: $sample_name\n" if $sample_name;
    print "Total number of mapped reads: $num_reads\n" if $num_reads;
    print "Total number of amplicons: ", scalar( keys %coverage_stats), "\n";
    print "Number of amplicons below the threshold: $low_total\n";
    printf "Percent of amplicons below the threshold: %.2f%%\n", ($low_total/scalar(keys %coverage_stats)) * 100;
    print "25%% Quartile Coverage: $quart1\n";
    print "50%% Quartile Coverage: $quart2\n";
    print "75%% Quartile Coverage: $quart3\n";
    
    print "total bases:          $total_bases\n";
    print "total non-zero bases: $total_nz_bases\n";
    print "total base reads:     $total_base_reads\n";
    print "mean base reads:      $mean_base_coverage\n";
    print "uniformity:           $uniformity\n";
    close $summary_fh;
}


sub get_metrics {
    my ($coverage_data,$total_base_reads) = @_;

    my ($total_nz_bases,$bases_over_mean);
    my $total_bases = keys %$coverage_data;
    my $mean_base_coverage = sprintf("%0.2f", $total_base_reads/$total_bases);
    for my $pos (keys %$coverage_data) {
        $total_nz_bases++ if $$coverage_data{$pos} != 0;
        $bases_over_mean++ if $$coverage_data{$pos} >= ($mean_base_coverage*0.2);
    }

    my $uniformity = sprintf("%0.2f", ($bases_over_mean/$total_bases)*100.00);

    #print "total bases:          $total_bases\n";
    #print "total non-zero bases: $total_nz_bases\n";
    #print "total base reads:     $total_base_reads\n";
    #print "mean base reads:      $mean_base_coverage\n";
    #print "uniformity:           $uniformity\n";
    return ($total_bases, $total_nz_bases, $total_base_reads, $mean_base_coverage, $uniformity);
}

sub gen_amplicon_coverage_tables {
    my $outdir = shift;
    my $lowcoverage_file = "$outdir/LowCoverageAmplicons.tsv";
    my $allcoverage_file = "$outdir/AllAmpliconsCoverage.tsv";
    open( my $aa_fh, ">", $allcoverage_file ) || die "Can't open the 'All Amplicons Coverage file for writing: $!";
    open( my $low_fh, ">", $lowcoverage_file ) || die "Can't open the Low Amplicons Coverage file for writing: $!";

    # Make tables of All Amplicon and Low Amplicon Coverage
    my $header = "Amplicon\tGene\tForward\tReverse\tFoward Proportion\tReverse Proportion\tMedian\tLength\n";
    print $aa_fh $header;
    print $low_fh $header;

    my $low_total = 0;
    for my $amplicon( keys %coverage_stats) {
        my ( $ampid, $gene, $length ) = split( /:/, $amplicon );
        my $outstring = join( "\t", $ampid, $gene, @{$coverage_stats{$amplicon}}, $length );
        print $aa_fh "$outstring\n";

        if ( $coverage_stats{$amplicon}[4] < $threshold ) {
            print $low_fh "$outstring\n";
            $low_total++;
        }
    }

    close $aa_fh;
    close $low_fh;
    return $low_total;
}

sub get_coverage_stats {
    my ($coverage_stats,$amp_coverage) = @_;
    for my $amplicon( sort keys %coverage_data ) {

        # amplicon strand info: length, median
        my $length = scalar(@{$coverage_data{$amplicon}->{'for'}});
        my @forward_reads = @{$coverage_data{$amplicon}->{'for'}};
        my @reverse_reads = @{$coverage_data{$amplicon}->{'rev'}};
        
        my $forward_median = median( \@forward_reads ); 
        my $reverse_median = median( \@reverse_reads );

        my @amp_coverage;
        for my $i ( 0..$#forward_reads ) {
            my $sum_fr = $forward_reads[$i] + $reverse_reads[$i]; 
            push( @all_coverage, $sum_fr );
            push( @amp_coverage, $sum_fr );
        }
        my $median = median( \@amp_coverage );

        my $total_reads = $forward_median + $reverse_median;

        my ( $forward_prop, $reverse_prop );
        if ( $total_reads != 0 ) {
            $forward_prop = sprintf( "%.3f", $forward_median/$total_reads );
            $reverse_prop = sprintf( "%.3f", 1 - $forward_prop );
        } else {
            $forward_prop = $reverse_prop = 0;
        }

        push( @{$coverage_stats{ join( ":", $amplicon, $length )}}, $forward_median, $reverse_median, $forward_prop, $reverse_prop, $median ); 
    }
    return;
}

sub get_base_coverage_data {
    my ($bambed, $regions_bed, $base_data) = @_;
    my $total_base_reads = 0;
    my $cmd = qq{ coveragebed -d -b $bambed -a $regionsbed };
    open(my $stream, "-|", $cmd);
    while (<$stream>) {
        my @fields = split;
        my $pos = "$fields[0]:" . ($fields[1] + ($fields[6]-1));
        $$base_data{$pos} = $fields[7];
    }
    return;
}

sub get_coverage_data {
    my ($bambed,$regionsbed,$coverage_data) = @_;

    # Use BEDtools to get amplicon coverage data for each amplicon
    my $get_forward_reads = qq{ grep "\\+\$" $bambed | coveragebed -d -b stdin -a $regionsbed };
    open( my $fcov, "-|", "$get_forward_reads" ) || die "Can't open the stream: $!";
    while (<$fcov>) {
        chomp;
        my @data = split;
        push( @{$coverage_data{join( ":", @data[3,5] )}->{'for'}}, $data[7] );
    }
    close $fcov;

    my $get_reverse_reads = qq{ grep "\\-\$" $bambed | coveragebed -d -b stdin -a $regionsbed };
    open( my $rcov, "-|", "$get_reverse_reads" ) || die "Can't open the stream: $!";
    while (<$rcov>) {
        chomp;
        my @data = split;
        push( @{$coverage_data{join( ":", @data[3,5] )}->{'rev'}}, $data[7] );
    }
    close $rcov;
    return;
}

sub median {
    # Create a median function in order to prevent having to use non-core modules
    my $data = shift;

    my @sorted_list = sort { $a <=> $b } @$data;
    my $length = @sorted_list;
    if ( $length%2 ) {
        return $sorted_list[int($length/2)];
    } else {
        return ($sorted_list[int($length/2)-1] + $sorted_list[int($length/2)])/2;
    }
}

sub quartile_coverage {
    # Homebrew Quartile function to avoid non-core modules
    my $data = shift;

    my @sorted_list = sort { $a <=> $b } @$data;

    my $second_quartile = median( $data );
    my @lower_half = grep { $_ < $second_quartile } @sorted_list;
    my @upper_half = grep { $_ > $second_quartile } @sorted_list;
    my $first_quartile = median( \@lower_half );
    my $third_quartile = median( \@upper_half );

    return( $first_quartile, $second_quartile, $third_quartile );    
}

sub proc_bed {
    my $input_bed = shift;
    (my $new_name = basename($$input_bed)) =~ s/(.*?)\.bed/$1_clean.bed/;
    my $output_bed = "$outdir/$new_name";

    open( my $fh, "<", $$input_bed ) || die "Can't open the BED file for processing: $!";
    open( my $out_fh, ">", $output_bed ) || die "Can't create the new BED file: $!";

    while (<$fh>) {
        if ( /track/ ) {
            print $out_fh $_;
            next;
        }

        my @fields = split;
        my ($gene, $pool);
        if ( $fields[4] eq '.' ) {
            ($gene, $pool) = $fields[-1] =~ /GENE_ID=(.*?);.*(Pool.?\d+)/;
        } else {
            $pool = $fields[4];
            $gene = $fields[5];
        }

        $pool =~ s/=//;

        print $out_fh join( "\t", @fields[0..3], $pool, $gene ), "\n";
    }
    close $fh;
    close $out_fh;

    return $output_bed;
}

sub gen_bam_bed {
    my $bamfile = shift;
    my $bambed = "$bamfile.bed";
    my $cmd = "bamToBed -i $bamfile";
    print "Generating a BED file from BAM $bamfile...";
    open( my $stream, "-|", $cmd );
    open( my $output_fh, ">", $bambed );

    while (<$stream>) {
        print {$output_fh} $_;
    }

    close $output_fh;
    print "Done!\n";
    return $bambed;
}
