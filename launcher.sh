#!/bin/bash
# Launcher script for AmpliconCoverageAnalysis-SA script
#
# TODO:
#    Need to write some code to clean up unecessary intermediate files:
#         Rplots.pdf
#         *clean.bed
#         BAM.bed?  => Should add cli arg to input bed file we already have one? So in that case keep it?
# 
# 8/8/2014 - D Sims
####################################################################################################################
VERSION="0.6.0_081014"
SCRIPTNAME=$(basename $0)
SCRIPTPATH=$(readlink -f $0)
SCRIPTDIR=$(dirname $SCRIPTPATH)

DEBUG=0

if [[ DEBUG -eq 1 ]]; then
    echo "#######  DEBUG  #######"
    echo "Script env info:"
    echo -e "\tscript name:   $SCRIPTNAME"
    echo -e "\tscript path:   $SCRIPTPATH"
    echo -e "\tscript dir:    $SCRIPTDIR"
    echo -e "#####################\n"
fi

USAGE="$(cat <<EOT
$SCRIPTNAME [options] <bamfile> <regions_bed> <sample_name>

Program to run the amplicon coverage analysis pipeline scripts.  
    -b    BAM BED file to use rather than making a new one
    -m    Minimum coverage threshold (DEFAULT: 450)
    -o    Output directory for all of the output data (DEFAULT: PWD)
    -v    Version information
    -h    Print this help text

EOT
)"

mincoverage=450  
outdir=$(pwd)
bambed=''

# Set up CLI opts
while getopts ":m:o:b:hv" OPT; do
    case "$OPT" in 
        b)
            bambed="$OPTARG"
            ;;
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

# Hi! My name is:
echo -e "\n::: Amplicon Coverage Analysis Pipeline :::\n"

# Check that we have all the args we need
if (( $# != 3 )); then
    echo "ERROR: Not enough arguments passed to script"
    echo "$USAGE"
    exit 1
else
    bamfile=$1
    if ! [[ -e "$bamfile" ]]; then
        echo "ERROR: The BAM file '$bamfile' does not exist"
        exit 1
    fi
    regions_bed=$2
    if ! [[ -e "$regions_bed" ]]; then
        echo "ERROR: The regions BED file '$regions_bed' does not exist"
        exit 1
    fi
    sample_name=$3
fi

if [[ DEBUG -eq 1 ]]; then
    echo "#######  DEBUG  #######"
    echo "CLI Opts:"
    echo -e "\tbam file:     $bamfile"
    echo -e "\tregions BED:  $regions_bed"
    echo -e "\tsample name:  $sample_name"
    echo -e "\tbam bed:      $bambed"
    echo -e "\tmin cov:      $mincoverage"
    echo -e "\toutdir:       $outdir"
    echo "#####################"
fi

run() {
    local exit_code=0
    
    #MSG=$( eval "$*" 2>&1 >/dev/null) 
    MSG=$( eval "$*" 2>&1 )
    exit_code=$?

    if [[ $exit_code != 0 ]]; then
        echo -e "ERROR: Nonzero exit '$exit_code' while running:\n \$ \`$*\`" >&2
        echo -e "\n$MSG"
        exit 1
    else
        echo "$MSG"
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
    test=(*clean.bed)
    declare -a temp_files=("Rplots.pdf" $outdir/*clean.bed "$outdir/$bambed")
    for file in "${temp_files[@]}"; do
        #echo rm -rf "$file"
        echo $(now) "Removing '$file'..."
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

# We're rolling...
echo "$(now) Starting pipeline..."

# Check the output directory
if ! [[ -d "$outdir" ]]; then
    echo -e "$(now) Directory '$outdir' does not exist.  Creating new directory\n"
    mkdir -p "$outdir"
fi
check_dir "$outdir"

# Set up and generate the necessary BAM BED file
if ! [[ $bambed ]]; then
    bambed=${bamfile/bam/bed}
    echo "$(now) Generating a BED file from '$bamfile'..."
    run "bamToBed -i $bamfile > \"${outdir}/$bambed\""
    echo "$(now) BED file '$bambed' generated successfully"
else
    echo "$(now) Using existing BED file '$bambed'..."
    if ! [[ -e $bambed ]]; then
        echo "ERROR: The BAM BED file '$bambed' does not exist"
        exit 1
    fi
fi

# Get the BAM size
echo  -e "\n$(now) Getting BAM file size..."
bamsize=$(run "samtools idxstats $bamfile | awk '{reads += \$3} END {print reads}'")
echo -e "\n$(now) $bamfile has $bamsize reads"

# Generate amp coverage tables
echo -e "\n$(now) Generating amplicon coverage tables..."
run "${SCRIPTDIR}/amplicon_coverage.pl -i -s $sample_name -t $mincoverage -r $bamsize -o $outdir $regions_bed $outdir/$bambed"
echo "$(now) Amplicon coverage tables successfully generated"

# Generate scatter plot
echo -e "\n$(now) Generating an amplicon coverage scatterplot..."
run "Rscript ${SCRIPTDIR}/coverage_scatter.R $sample_name $mincoverage $outdir"
echo "$(now) Scatter plots generated successfully"

# Generate strand bias tables
echo -e "\n$(now) Generating amplicon bias plots..."
run "Rscript ${SCRIPTDIR}/strand_coverage.R $outdir/AllAmpliconsCoverage.tsv $outdir"
echo "$(now) Bias plots successfully generated"

# Clean up
echo -e "\n$(now) Cleaning up intermediate files..."
cleanup
echo -e "\n$(now) Done!"
mv $logfile $outdir
