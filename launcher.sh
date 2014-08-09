#!/bin/bash
# Launcher script for AmpliconCoverageAnalysis-SA script
#
# TODO:
#    Need to set up package environment so that we can find all the scripts
#    Need to write some code to clean up unecessary intermediate files:
#         Rplots.pdf
#         *clean.bed
#         BAM.bed?  => Should add cli arg to input bed file we already have one? So in that case keep it?
# 
# 8/8/2014 - D Sims
####################################################################################################################
VERSION="0.4.0_080914"
SCRIPTNAME=$(basename $0)
SCRIPTPATH=$(readlink -f $0)
SCRIPTDIR=$(dirname $SCRIPTPATH)

DEBUG=1

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
    -m    Minimum coverage threshold (DEFAULT: 450)
    -o    Output directory for all of the output data (DEFAULT: PWD)
    -v    Version information
    -h    Print this help text

EOT
)"

declare -a required_progs=("bamToBed" "samtools") 

mincoverage=450  
outdir=$(pwd)

# Set up CLI opts
while getopts m:o:hv OPT; do
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
            help
            exit 1
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            help
            exit 1
            ;;
    esac
done
shift $((OPTIND - 1 ))

# Check that we have all the args we need
if (( $# != 3 )); then
    echo "ERROR: Not enough arguments passed to script"
    echo "$USAGE"
    exit 1
else
    bamfile=$1
    if ! [[ -e "$bamfile" ]]; then
        echo "ERROR: BAM file '$bamfile' does not exist"
        exit 1
    fi
    regions_bed=$2
    if ! [[ -e "$regions_bed" ]]; then
        echo "ERROR: Regions BED file '$regions_bed' does not exist"
        exit 1
    fi
    sample_name=$3
fi

bambed=${bamfile/bam/bed}

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
    continue
    declare -a temp_files=("Rplots.pdf" "*clean.bed" "$bambed")
    for file in $temp_files; do
        echo rm -rf "$file"
    done
}

check_dir() {
    shopt -s nullglob
    outdir=$1
    re='[a-zA-Z0-9]+\.(tsv|txt|bed)$'

    dir_contents=("$outdir"/*)

    if (( "${#dir_contents}" )); then
        for file in $dir_contents; do
            echo "$file"
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

check_env ${required_progs[@]}

# Set up a log file
logfile="$(date +"%Y%m%d").aac.log.txt"
now() { now=$(date +%c); echo -n "[$now]:";}
exec > >(tee $logfile)
exec 2>&1

# Check the output directory
if ! [[ -d "$outdir" ]]; then
    echo -n "$(now) Directory '$outdir' does not exist.  Creating new directory\n\n"
    mkdir -p "$outdir"
fi
check_dir "$outdir"
exit 0

# Generate a BAM.bed file.
echo "$(now) Generating a BED file from '$bamfile'..."
run "bamToBed -i $bamfile > \"${outdir}/$bambed\""
echo "$(now) BED file '$bambed' generated successfully"

# Get the BAM size
echo  "\n$(now) Getting BAM file size..."
bamsize=$(run "samtools idxstats $bamfile | awk '{reads += \$3} END {print reads}'")
echo "$(now) $bamfile has $bamsize reads"

# Generate amp coverage tables
echo "\n$(now) Generating amplicon coverage tables..."
run "${SCRIPTDIR}/amplicon_coverage.pl -i -s $sample_name -t $mincoverage -r $bamsize -o $outdir $regions_bed $bambed"
echo "$(now) Amplicon coverage tables successfully generated"

# Generate scatter plot
echo "\n$(now) Generating an amplicon coverage scatterplot..."
run "Rscript ${SCRIPTDIR}/coverage_scatter.R $sample_name $mincoverage $outdir"
echo "$(now) Scatter plots generated successfully"

# Generate strand bias tables
echo "\n$(now) Generating amplicon bias plots..."
run "Rscript ${SCRIPTDIR}/strand_coverage.R $outdir/AllAmpliconsCoverage.tsv $mincoverage $outdir"
echo "$(now) Bias plots successfully generated"

# Clean up
echo "\n$(now) Cleaning up intermediate files..."
cleanup
echo "\n$(now) Done!"
mv $logfile $outdir


