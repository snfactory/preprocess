#!/bin/sh
######################################################################
## Filename:      quick_preprocess.sh
## Version:       $Revision$    
## Author:        $Author$
## Id:            $Id$
######################################################################

scriptname=$(basename $0)
cvsname=$(echo '$Name$' | sed 's/.*: //;s/..$//')
version=$(echo '$Revision$' | sed 's/.*: //;s/..$//')
version=$(echo "${cvsname:=developer}-$version")

# Default reference files
biasm=biasmodel
bias=bias
darkm=darkmodel
dark=dark

# ============================================================
# Usage
function usage() {
    cat <<EOF
Program: quick_preprocess Version $version

Available options:
[-c channel] -i incat [-o outcat] [-b biasmodel] [-B biasmap]
   [-d darkmodel] [-D darkmap]  [-t timeon ] [-T timeonfile] [-C] [-F]
   [-f] [-D] [-h] [-v]

where:
    channel : Channel (R or B).
              Default: automatic detection.
      incat : Input frame or catalog. If the catalog is too long, it will be splitted
              into smaller catalogs
    outcat  : Output frame or catalog default is to propduce a P file
  biasmodel : use the specified bias model
       bias : Use specified bias map
  darkmodel : use the specified dark model
       dark : use the specified dark map
     timeon : the time since on value to be used for a single file
              for multi-file processing, use timeonfile
 timeonfile : list of the times on. The format is 'file time' per line 
              with no header
              look in the preproc function if you want to get it right
     C-flag : local copy of input files. Useful if input files are write 
              protected.
     F-flag : full processing of P channel.
     f-flag : no confirmation asked
     D-flag : debug mode
     h-flag : displays help and exits
     v-flag : displays version and exits

Purpose:
Long preprocessing, with all substractions possible

Remarks:
- The automatic detection of the input channel is still very dumb: if inframe
  is of the form *_X.fits , it corresponds to X channel.
- If the input catalog is too long, it will be splitted in smaller parts
- The file names paths shall not contain a '.' (this is a general IFU 
  limitation)
- WARNING: the TIMEON qualifier has to be stored by someone inside the 
  exposure.
  If the T or t option are provided, it is used
  Else if the qualifier exists, it is assumed to be correct. 
  Else the DB will be used to retreive it.
- The time on file may be longer than the actual catalog. it just has to 
  contain the file
EOF
echo '$Id$'
}

# ============================================================
# Fatal error termination in the script 
# Usage: die status $LINENO "message"
function die() {
    echo "ERROR $1 in $scriptname (l.$2): $3"
    exit $1
}

# ============================================================
# Preprocess core routine
# Usage: preproc
# it relies heavily on environment variables set in main
function preproc() {

    #CP the input files if required
    if [ $cpinput ] ; then
	for file in `sed 1d $tmptmpincat` ; do
	    cp $file ./
	done
	sed 's%.*/%%' $tmptmpincat > tmp3in.cat
    else
	cp $tmptmpincat tmp3in.cat
    fi

    # Plug the timeon if needed
    for file in $(sed "1d" tmp3in.cat ); do
	timedesc=`rd_desc -in $file -desc TIMEON -quiet`
	darktime=`rd_desc -in $file -desc DARKTIME -quiet`
	[ $DEBUG ] && echo "file $file timeon file $timeon time for this file $timedesc darktime $darktime"
	if [ $timeon ]  || [ $timeonval ]; then
	    if [ -z $timeonval ] ; then
		timeonval=`grep $file $timeon | sed 's/.* //'`
	    fi
	    if [ $DEBUG ] ; then
		echo wr_desc -in $file -desc TIMEON -val $(echo $darktime + $timeonval | bc ) -quiet
	    fi
	    wr_desc -in $file -desc TIMEON -val $(echo $darktime + $timeonval | bc ) -quiet || die 1 $LINENO "not possible to write the timeon descriptor to $file"
	elif [ -z "$timedesc" ] ; then
	    timeonval=`try_db33 $file | sed 's/.* //'`
	    if [ $DEBUG ] ; then
		echo "wr_desc -in $file -desc TIMEON -val $(echo $darktime + $timeonval | bc ) -quiet" 
	    fi
	    wr_desc -in $file -desc TIMEON -val $(echo $darktime + $timeonval | bc ) -quiet || die 1 $LINENO "not possible to write the timeon descriptor to $file"
	else 
	    [ $DEBUG ] && echo "already a time on $timedesc nothing done"
	fi
	timeonval=""
    done

    # build the output catalog if needed
    if [ -z $outcat ] ; then
	sed 's%.*/%% ; /.fits/s/^/P/ ' tmp3in.cat > $tmptmpoutcat
    fi

    if [ $DEBUG ] ; then 
	echo preprocess -in tmp3in.cat -out $tmptmpoutcat $option $quiet $noask $full

    fi
    preprocess -in tmp3in.cat -out $tmptmpoutcat $option $quiet $noask $full
}

# ============================================================
# Check if a frame is a raster, and embed it in full-frame if needed.
# Usage: check_and_embed_raster file || echo "Un-embedable"
function is_raster() {
    local raster
    # Check if it's a raster
    raster=$(echo "$(stat -L -c %s $1) < 1.6*10^7" | bc) # Less than 16Mb
    echo $raster
}


# ############################################################

# Options ==============================

# Parser ..............................
while getopts "c:i:o:b:B:d:D:T:t:CFfgvh" OPTION ; do
    case "$OPTION" in
        c) channel="$OPTARG" ;; # Input channel
        i) incat="$OPTARG" ;; # Input frame
        o) outcat="$OPTARG" ;; # Output frame
        b) biasmodel="$OPTARG" ;;  # Bias model
        B) biasmap="$OPTARG" ;;    # Bias map
	d) darkmodel="$OPTARG" ;;    # Dark model
        D) darkmap="$OPTARG" ;; # Dark map
        t) timeonval="$OPTARG" ;; # The unique time on
        T) timeon="$OPTARG" ;; # The list of time on
        C) cpinput="yes" ;;    # Cp input files in the cwd before processing
        f) noask="-noask" ;;    # Do not ask before overwriting files
        F) full="-all" ;;    # Process the 4 amps of P channel
        g) DEBUG="-d"; echo "*** DEBUG mode ***" ;; 
        v) echo "$version" ; exit 0 ;;
        h) usage; exit 0 ;;
        ?) usage; exit 1 ;;
    esac
done
[ $DEBUG ] || quiet="-quiet"

echo "===== $(basename $0) $version -- $(date) ====="

if [ -z $incat ] ; then
    die 1 $LINENO "Arguments of options -i is missing"
fi
if [ ! -r $incat ] ; then
    die 1 $LINENO "Cannot read input files $incat "
fi
if [ $outcat ] && [ ! -r $outcat ]; then
    die 1 $LINENO "Cannot read output file $outcat"
fi
if [ $biasmap ] && [ ! -r $biasmap ]; then
    die 1 $LINENO "Cannot read bias file $bias"
fi
if [ $biasmodel ] && [ ! -r $biasmodel ]; then
    die 1 $LINENO "Cannot read bias model $biasmodel"
fi
if [ $darkmap ] && [ ! -r $darkmap ]; then
    die 1 $LINENO "Cannot read bias file $bias"
fi
if [ $darkmodel ] && [ ! -r $darkmodel ]; then
    die 1 $LINENO "Cannot read dark model $darkmodel"
fi
if [ $timeon ] && [ $timeonval ]; then
    die 1 $LINENO "Cannot use both -t and -T options"
fi
if [ $timeon ] && [ ! -r $timeon ]; then
    die 1 $LINENO "Cannot read timeon file $timeon"
fi


datadir="$(dirname $0)/../data"
[ -d $datadir ] || \
    die 1 $LINENO "Cannot find reference calibration directory $datadir"

# check length of catalogs : very basic matching .fits is a file .cat is a catalog
case $incat in
    *.fits) ;;
    *.cat) iscat=yes ;;
    *) die 1 $LINENO "Input file doesn't look like a known format. use .cat of .file" ;;
esac

# outcat is expected to be of the same type as incat
if [ $iscat ] ; then
    length=$(wc -l $incat | sed 's/ .*//')
    length=$(( $length - 1 ))
else
    length=1
fi

if [  "$iscat" -a "$outcat" ] ; then
    if [ $(wc -l $incat | sed 's/ .*//') -ne $(wc -l $outcat | sed 's/ .*//') ] ; then
	die 1 $LINENO "input catalog and output catalog are not of the same length"
    fi
fi

# Automatic detection of the input channel ..............................
if [ -z $channel ]; then
    if [ $iscat ] ; then
	inframe=`tail -n 1 $incat`
    else
	inframe=$incat
    fi
    # Should be done with descriptors, for the moment with filename
    case $inframe in
        *_B.fits) channel=B ;;
        *_R.fits) channel=R ;;
        *_P.fits) channel=P ;;
        *) die 1 $LINENO "Automatic channel detection failed (see option -c)" ;;
    esac
fi
echo "Input channel: $channel"


# in a catalog, all files should be of the same channel and raster
israster=`is_raster $inframe`

if [ $iscat ] ; then
    for file in `sed "1d" $incat` ; do
	if [ $(expr "$file" : ".*_${channel}.fits" ) -eq 0 ] ; then
	    die 1 $LINENO "file $file is not of channel $channel"
	fi
	if [ `is_raster $file` != $israster ] ; then
	    die 1 $LINENO "file $file is not of same raster as $inframe"
	fi
    done
fi

# Default values  ==============================
[ $DEBUG ] && echo channel is $channel raster is $israster frame is $inframe
if [ $israster -eq 0  ] && [ $channel != "P" ] ; then
    if [ "( -z $biasmodel ) -a ( $israster -eq 0 )" ] ; then
	biasmodel=${biasm}_${channel}.txt
	ln -s $datadir/$biasmodel ./
    fi
    if [ "( -z $biasmap ) -a ( $israster -eq 0 )" ] ; then
	biasmap=${bias}_${channel}.fits
	ln -s $datadir/$biasmap ./
    fi
    if [ "(-z $darkmodel ) -a ( $israster -eq 0 ) " ] ; then
	darkmodel=${darkm}_${channel}.txt
	ln -s $datadir/$darkmodel ./
    fi
    if [ "( -z $darkmap ) -a ( $israster -eq 0 )" ] ; then
	darkmap=${dark}_${channel}.fits
	ln -s $datadir/$darkmap ./
    fi
    if [ "( $israster -eq 0 ) -a ( $channel != "P" )" ] ; then
	option=" -bm $biasmodel -bias $biasmap -dm $darkmodel -dark $darkmap"
    fi
fi
echo "option is $option"

# Catalog splitting (if needed) ==============================

tmpincat=tmpin.cat
tmptmpincat=tmptmpin.cat
tmptmpoutcat=tmptmpout.cat

if [ $outcat ] ; then
    tmpoutcat=tmpout.cat
    cp $outcat $tmpoutcat
fi

if [ -z $iscat ] ; then
    cat <<EOF > $tmpincat
######## a temp catalog
EOF
    echo $incat >> $tmpincat
else 
    cp $incat $tmpincat
fi

# preprocessing the catalog
while [ `wc -l $tmpincat | sed "s/ .*//"` -gt 20 ] ; do
    head -n 21 $tmpincat > $tmptmpincat
    if [ $outcat ] ; then
	head -n 21 $tmpoutcat > $tmptmpoutcat
    fi
    preproc 
    sed "2,21d" $tmpincat > $tmptmpincat
    mv $tmptmpincat $tmpincat
    if [ $outcat ] ; then
	sed "2,21d" $tmpoutcat > $tmptmpoutcat
	mv $tmptmpoutcat $tmpoutcat
    fi
done

cp $tmpincat $tmptmpincat
if [ $outcat ] ; then
    cp $tmpoutcap  $tmptmpoutcat
fi
preproc 

[ $DEBUG ] || rm -f tmp*

exit 0
