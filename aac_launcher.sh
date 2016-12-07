#!/bin/bash
# Launcher script for AmpliconCoverageAnalysis-SA script
#
# TODO:
#    - Add some better error checking.
# 
# 8/8/2014 - D Sims
####################################################################################################################
VERSION="1.5.4_120716"
SCRIPTNAME=$(basename $0)
SCRIPTDIR=$(dirname $(readlink -f $0))

USAGE="$(cat <<EOT
$SCRIPTNAME [options] <bamfile> <regions_bed> <sample_name>

Program to run the amplicon coverage analysis pipeline scripts on a BAM file from the Ion Torrent platform. Required
input is a BAM file, a regions BED file, and a sample name for formtting the output. Note that it is important to have
a properly formatted BED file, or else the pipeline will not run as the necessary data will not be extracted.

    -m    Minimum coverage threshold (DEFAULT: 450)
    -o    Output directory for all of the output data (DEFAULT: PWD)
    -v    Version information
    -h    Print this help text

EOT
)"

mincoverage=450  
outdir=$(pwd)

# Set up CLI opts
while getopts ":m:o:hv" OPT; do
    case "$OPT" in 
        m)
            if ! [[ $OPTARG =~ ^[0-9]+$ ]]; then
                echo "ERROR: '$OPTARG' not valid for opt '-m'. Minimum coverage must be a number"
                exit 1
            else
                mincoverage=$OPTARG
            fi
            ;;
        o)
            outdir=$OPTARG
            ;;
        h)
            echo -e "$SCRIPTNAME - $VERSION\n"
            echo "$USAGE"
            exit 0
            ;;
        v)
            echo "$SCRIPTNAME - $VERSION"
            exit 0
            ;;
        :)
            echo "Option "-$OPTARG" requires an argument\n"
            echo "$USAGE"
            exit 1
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            echo "$USAGE"
            exit 1
            ;;
    esac
done
shift $((OPTIND - 1 ))

filecheck() {
    for file in $@; do 
        if ! [[ -e $file ]]; then 
            echo "ERROR: the file '$file' does not exist!"
            echo
            exit 1
        fi
    done
}

run() {
    local exit_code=0
    
    MSG=$( eval "$*" 2>&1 )
    exit_code=$?

    if [[ $exit_code != 0 ]]; then
        echo -e "\nERROR: Nonzero exit '$exit_code' while running:\n \$ \`$*\`" >&2
        echo -e "\n$MSG"
        exit 1
    else
        echo -n "$MSG"
    fi
}

check_env() {
    # Check to be sure we have the required external progs
    required_progs=("$@")

    for prog in "${required_progs[@]}"; do
        if ! command -v $prog >/dev/null 2>&1; then
            echo "ERROR: I can't find the required program '$prog'. Is it installed?"
            exit 1;
        fi
    done
}

cleanup() {
    declare -a temp_files=("Rplots.pdf" $outdir/*clean.bed "$outdir/$bambed")
    for file in "${temp_files[@]}"; do
        echo -e "\tRemoving '$file'..."
        rm -rf "$file"
    done
}

check_dir() {
    shopt -s nullglob
    outdir=$1
    re='[a-zA-Z0-9]+\.(tsv|txt|bed)$'

    dir_contents=("$outdir"/*)

    if (( "${#dir_contents}" )); then
        for file in $dir_contents; do
            if [[ $file =~ $re ]]; then
                echo "WARNING: Output directory '$outdir' not empty and may have data that will be overwritten!"
                echo -n "(c)ontinue (e)xit (n)ew directory? "
                while [[ 1 ]]; do
                    read resp
                    case "$resp" in 
                        c)
                            echo "Continuing..."
                            return 1 
                            ;;
                        e)
                            echo "Exiting!"
                            exit 1
                            ;;
                        n)
                            echo -n "New dir: "
                            read dir
                            mkdir -p "$dir"
                            check_dir "$dir"
                            return 1
                            ;;
                        *)
                            echo "ERROR: '$resp' not a valid option"
                            echo -n "(c)ontinue (e)xit (n)ew directory? "
                            ;;
                    esac
                done
            fi
        done
    fi
}

# Check that we have the required helper progs to run this package
declare -a required_progs=("bamToBed" "samtools") 
check_env ${required_progs[@]}

# Set up a log file
logfile="$(date +"%Y%m%d").aac.log.txt"
now() { now=$(date +"%m-%d-%Y %T"); echo -n "[$now]:";}
exec > >(tee $logfile)
exec 2>&1

# Hi! My name is:
echo -e "\n::: Amplicon Coverage Analysis Pipeline :::\n"

# Check that we have all the args we need
if (( $# < 3 )); then
    echo "ERROR: Not enough arguments passed to script"
    echo "$USAGE"
    exit 1
else
    bamfile=$1
    regions_bed=$2
    sample_name=$3
    filecheck $bamfile $regions_bed
fi

# We're rolling...
echo "$(now) Starting pipeline $bamfile..."
echo "$(now) Params as passed to the pipeline:"
echo -e "\tMinimum coverage  => $mincoverage"
echo -e "\tBAM BED File      => $bambed"
echo -e "\tOutput Directory  => $outdir"
echo -e "\tRegions BED File  => $regions_bed"
echo -e "\tBAM File          => $bamfile"
echo -e "\tSample name       => $sample_name"

# Check the output directory
if ! [[ -d "$outdir" ]]; then
    echo "$(now) Directory '$outdir' does not exist.  Creating new directory."
    mkdir -p "$outdir"
fi
check_dir "$outdir"

# Set up and generate the necessary BAM BED file
bambed="${bamfile}.bed"
echo -n "$(now) Generating a BED file from '$bamfile'..."
run "bamToBed -i $bamfile > \"${outdir}/$bambed\""
echo "Done!"
echo "$(now) BED file '$bambed' generated successfully"

# Get the BAM size
echo "$(now) Getting BAM file size..."
if ! [[ -e ${bamfile}.bai ]]; then 
    echo -n "$(now)     BAM file '$bamfile' is not indexed!  Indexing the file now..."
    run "samtools index $bamfile"
    echo "Done!"
fi
bamsize=$(run "samtools idxstats $bamfile | awk '{reads += \$3} END {print reads}'")
echo -e "$(now)     $bamfile has $bamsize reads"

# Generate amp coverage tables
echo -n "$(now) Generating amplicon coverage tables..."
run "${SCRIPTDIR}/amplicon_coverage.pl -i -s $sample_name -t $mincoverage -r $bamsize -o $outdir $regions_bed -b $outdir/$bambed 2>&1 > /dev/null"
echo "Done!"
echo "$(now) Amplicon coverage tables successfully generated"

# Generate scatter plot
echo -n "$(now) Generating an amplicon coverage scatterplot..."
run "Rscript ${SCRIPTDIR}/coverage_scatter.R $sample_name $mincoverage $outdir 2>&1 > /dev/null"
echo "Done!"
echo "$(now) Scatter plots generated successfully"

# Generate strand bias tables
echo "$(now) Generating amplicon bias plots..."
run "Rscript ${SCRIPTDIR}/strand_coverage.R $outdir/AllAmpliconsCoverage.tsv $outdir 2>&1 > /dev/null"
echo "Done!"
echo "$(now) Bias plots successfully generated"

# Clean up
echo "$(now) Cleaning up intermediate files..."
cleanup
echo "$(now) Done!"
mv $logfile $outdir
echo "Pipeline complete.  Data ready to view in directory: $outdir"
