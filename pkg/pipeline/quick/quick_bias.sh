#!/bin/sh
######################################################################
## Filename:      quick_bias.sh
## Version:       $Revision$    
## Author:        $Author$
## Id:            $Id$
######################################################################

# arg : we should have a list of files
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
 -o outbias [-m maximages] [-f][-d][-h][-v] bias1.fits bias2.fits [bias*.fits]

where:
    outbias : Output bias name.
  maximages : The maximum number of images which may be simultaneously
              present in the memory.
	      The higher, the better performance - but the higher chance to
	      overflow the RAM.
	      Note that it is not the max number of input images. This one can 
              be very high (up to 1000 should go fine).
	      default : 7 (corresponding to approx 1GB RAM)
     f-flag : Do not ask before overwriting files.
     d-flag : Debug mode

Example:

    quick_bias -o bias_B.fits -f *_24_B.fits

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
while getopts "o:m:fdvh" OPTION ; do
    case "$OPTION" in
        o) outbias="$OPTARG" ;; # Output datacube
	m) maximages="$OPTARG" ;; # max space for processing
        f) noask="-noask" ;;    # Do not ask before overwriting files
	d) DEBUG="-d";echo "*** DEBUG mode ***";;
        v) echo "$version" ; exit 0 ;;
        h) usage; exit 0 ;;
        ?) usage; exit 1 ;;
    esac
done

shift $(($OPTIND-1))
inframes=$@

if [ -z "$inframes" ] || [ -z $outbias ]; then
    die 1 $LINENO "Arguments of options -o missing or no input files"
fi
if [ -e Pfiles.cat ] ; then
    if [ -n "$noask" ] ; then
	rm -f Pfiles.cat
    else 
	die 1 $LINENO "Pfiles.cat already exists - is a concurrent processing runing on this directory?"
    fi
fi
### Translate the number of images in number of lines for 1 image.
# We have to hold in fact ALL images in memory.
nlines=$(echo "4096*$maximages/($(echo $inframes | wc -w)+1)" | bc)

cat <<EOF > Pfiles.cat
#### Processed files catalog
EOF

### Actual processing here

for file in $inframes ; do
    Pfile=`echo $file | sed 's%.*/%% ; s%^%P%'`
    echo $Pfile >> Pfiles.cat
    [ $DEBUG ] && echo "preprocess -in $file -out $Pfile"
    preprocess -in $file -out $Pfile $noask
done
[ $DEBUG ] && echo "stack_image -in Pfiles.cat -out $outbias -nlines $nlines"
stack_image -in Pfiles.cat -out $outbias -nlines $nlines $noask
rm -f Pfiles.cat

exit 0
