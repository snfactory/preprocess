#!/bin/sh
######################################################################
## Filename:      quick_bias.sh
## Version:       $Revision$    
## Author:        $Author$
## Id:            $Id$
######################################################################

# arg : we should have a list of files (or a catalog)
# they will be processed
# then a bias will be computed out of them

version=$(echo '$Revision$' | sed 's/.*: //;s/..$//')
maximages=7

# ===================================================
# Usage
function usage() {
    cat <<EOF
Program: quick_bias  Version $version

Available options:
 -o outbias [-C incat] [-m maximages] [-f][-d][-h][-v] [bias1.fits [bias2.fits ...]]

where:
    outbias : Output bias name.
      incat : Input catalog. Can replace multiple arguments.
  maximages : The maximum number of images which may be simultaneously present
              in the memory.  The higher, the better performance - but the
              higher chance to overflow the RAM.  Note that it is not the max
              number of input images. This one can be very high (up to 1000
              should go fine).
	      default : 7 (corresponding to approx 1GB RAM)
     f-flag : Do not ask before overwriting files.
     d-flag : Debug mode

Examples:
    quick_bias -o bias_B.fits -f *_24_B.fits
    quick_bias -o bias_B.fits -C mybiases_B.cat
EOF
}

# ============================================================
# Fatal error termination in the script 
# Usage: die status $LINENO "message"
function die() {
    echo "ERROR $1 in $scriptname (l.$2): $3"
    exit $1
}


# ############################################################

# Options ==============================

# Parser ..............................
while getopts "o:m:C:fdvh" OPTION ; do
    case "$OPTION" in
        o) outbias="$OPTARG" ;; # Output datacube
	m) maximages="$OPTARG" ;; # max space for processing
	C) catalog="$OPTARG" ;;   # Input catalog
        f) noask="-noask" ;;    # Do not ask before overwriting files
	d) DEBUG="-d";echo "*** DEBUG mode ***";;
        v) echo "$version" ; exit 0 ;;
        h) usage; exit 0 ;;
        ?) usage; exit 1 ;;
    esac
done

# Read arguments
shift $(($OPTIND-1))
inframes=$@

if [ -z "$inframes" ] && [ -z $catalog ]; then
    die 1 $LINENO "No input files or option -C missing."
fi
if [ -z $outbias ]; then
    die 1 $LINENO "Option -o missing."
fi

# Create input catalog if needed
if [ -z "$catalog" ]; then
    incat=qb_$$.cat
    echo "# quick_bias $version -- $(date)" > $incat
    echo $inframes | tr ' ' '\n' >> $incat
    [ $DEBUG ] && echo "Creating input catalog $incat"
else
    incat=$catalog
    [ $DEBUG ] && echo "Using input catalog $incat"
fi

### Translate the number of images in number of lines for 1 image.
# We have to hold in fact ALL images in memory.
nframes=$(cat $incat | grep -v '^#' | wc -l)
nlines=$(echo "4096*$maximages/($nframes+1)" | bc)

# Create preprocess output catalog
outcat=Pqb_$$.cat
[ $DEBUG ] && echo "Creating output catalog $outcat"
cat $incat | sed -e '1! s%.*/%P%' > $outcat

### Actual processing here

[ $DEBUG ] && echo "preprocess -in $incat -out $outcat"
for infile in `sed '1 d' $incat`; do
    outfile=`echo $infile | sed -e '1! s%.*/%P%'`
    preprocess -in $infile -out $outfile $noask
done

[ $DEBUG ] && echo "stack_image -in $outcat -out $outbias -nlines $nlines"
stack_image -in $outcat -out $outbias -nlines $nlines $noask

# Cleaning
if [ $DEBUG ]; then
    [ -z "$catalog" ] && echo "DEBUG: Catalog $incat not removed."
    echo "DEBUG: Catalog $outcat not removed."
else
    [ -z "$catalog" ] && rm -f $incat
    rm -f $outcat
fi

exit 0
