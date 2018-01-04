#!/bin/bash

# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#  
# A copy of the GNU General Public License is included in the file COPYING and is
# available at http://www.gnu.org/licenses/.
#
# Copyright 2015-2017 Hyrum D. Carroll

DEBUG=1

NUM_TOTAL_ITERATIONS=5
numTotalIterations=$NUM_TOTAL_ITERATIONS

if [ -z "$1" ]; then
    echo "Usage: $0 <directory of directories that have retrieval results> <timing output file base name>" &>2
    exit 1;
fi

runsDir=$1
shift
if [ $DEBUG -eq 1 ]; then echo "DEBUG: runsDir: $runsDir"; fi

timingOutBase="$1"
shift
if [ $DEBUG -eq 1 ]; then echo "DEBUG: timingOutBase: $timingOutBase"; fi


timingFileDir="${runsDir%/}"
timingFileDir="${timingFileDir%/*}"

cmd="cd $timingFileDir"
echo $cmd
$cmd

ks="1 3 5 20"
allTimingOutFile="$timingOutBase-allTiming.csv"
echo "allTimingOutFile: $allTimingOutFile"
timingOutFile="$timingOutBase.timing.out"
nrIterTimingsStr="3 k 0"
finalTimingsStr="3 k 0"

if [ ! -s "$allTimingOutFile" ]; then
    if [ $DEBUG -eq 1 ]; then echo "DEBUG: Gathering timing info . . ."; fi
    headerStr="#tag,dataset"
    iter=1
    while [ $iter -lt $NUM_TOTAL_ITERATIONS ]; do
	headerStr="$headerStr,iter$iter"
	iter=$(( iter + 1 ))
    done
    headerStr="$headerStr,iterLast"
    headerStr="$headerStr,finalTiming,sumTiming"
    echo "$headerStr" > $allTimingOutFile
    for dir in `basename $runsDir`/*; do
	if [ ! -d "$dir" ]; then
	    continue
	fi

	# collect iterativeDBTiming (real) timing from $app.<num>.nr.out
	# collect final (real) timing from $app.<num>.final.out
	# sum timings
	
	nrIterTiming=""  # all of the (real) timing for each iteration, separated by commas
	nrTimings=""     # all of the (real) timing for each iteration, formated for bc
	finalTiming=""
	sumTiming=""

	baseDir=`basename $dir`
	dataset="${baseDir%.*}"
	tag="${baseDir#*.}"
	app="${tag%%-*}"
	# if [ $DEBUG -eq 1 ]; then echo "tag: $tag; app: $app" >&2; fi
	if [ "$app" == "hmmsearch" ]; then
	    app="hmmer"
	elif [ "$tag" == "psiblast-non_iterative" ]; then
	    app="$tag"
	elif [[ "$tag" == psisemiglobal-non_iterative* ]]; then
	    app="psisemiglobal-non_iterative"
	    if [ $DEBUG -eq 1 ]; then echo "app: $app (psisemiglobal-non_iterative)" >&2; fi
	elif [[ "$tag" == psisemiglobal* ]]; then
	    # PSI-SemiGLOBAL can save the PSSM and then exit, requiring one less iteration
	    numTotalIterations=$(( NUM_TOTAL_ITERATIONS - 1 ))
	fi
	
	iter=1
	while [ $iter -le $NUM_TOTAL_ITERATIONS ]; do
	    nrTiming=""
	    # insert blank entries for iterations that were skipped (i.e. for PSI-SemiGLOBAL)
	    if [ $numTotalIterations -ne $NUM_TOTAL_ITERATIONS ]; then
		if [ $iter -eq $numTotalIterations ]; then
		    while [ $iter -lt $NUM_TOTAL_ITERATIONS ]; do
			nrIterTiming="$nrIterTiming,$nrTiming"
			iter=$(( iter + 1 ))
		    done
		fi
	    fi
	    nrIterOutputFile="$dir/$dataset.$app.iterationsDb.out"
	    # the last iteration is special and does not have the iteration number appended to it
	    if [ $iter -lt $numTotalIterations ]; then
		nrIterOutputFile="$nrIterOutputFile$iter"
	    fi
	    if [ -f "$nrIterOutputFile" ]; then
		gzipped=0
	    elif [ -f "$nrIterOutputFile.gz" ]; then
		gzipped=1
		#echo "gunzipping $nrIterOutputFile.gz"
		gunzip "$nrIterOutputFile.gz"
		if [ ! -f "$nrIterOutputFile" ]; then
		    echo "ERROR: nrIterOutputFile, \"$nrIterOutputFile\" . . ." >&2
		    ls -la "$PWD/$nrIterOutputFile" >&2
		    exit 1
		fi
	    fi
	    # echo "nrIterOutputFile: $nrIterOutputFile (iter: $iter)"
	    if [ -s "$nrIterOutputFile" ]; then
		nrTiming=`grep "real [0-9.]" $nrIterOutputFile`
		nrTiming="${nrTiming##*real }"
		# echo "nrTiming: $nrTiming"
		if [ -n "$nrTiming" ]; then
		    nrTimings="$nrTimings $nrTiming + "
		    nrIterTimingsStr="$nrIterTimingsStr $nrTiming + "
		else
		    echo "ALERT: Did not find timing information in $nrIterOutputFile! (gzipped: $gzipped)" >&2
		fi
		#echo ",,,$nrTiming"

		if [ $gzipped -eq 1 ]; then
		    rm -f "$nrIterOutputFile.gz"
		    gzip $nrIterOutputFile
		fi
	    fi
	    nrIterTiming="$nrIterTiming,$nrTiming"
	    # echo "nrIterTiming: $nrIterTiming"
	    iter=$(( iter + 1 ))
	done

	finalOutputFile="$dir/$dataset.$app.final.out"
	# echo "finalOutputFile: $finalOutputFile"
	gzipped=0
	if [ -f "$finalOutputFile.gz" ]; then
	    gzipped=1
	    #echo "gunzipping $finalOutputFile.gz"
	    gunzip "$finalOutputFile.gz"
	    if [ ! -f "$finalOutputFile" ]; then
		echo "ERROR: finalOutputFile, \"$finalOutputFile\" . . ." >&2
		ls -la "$PWD/$finalOutputFile" >&2
		exit 1
	    fi
	fi
	if [ -s "$finalOutputFile" ]; then
	    finalTiming=`grep "real [0-9.]" $finalOutputFile`
	    #echo "finalTiming:,,,$nrIterTiming,$finalTiming,$sumTiming"
	    finalTiming="${finalTiming##*real }"
	    #echo "finalTiming:,,,$nrIterTiming,$finalTiming,$sumTiming"
	    if [ -n "$finalTiming" ]; then 
		finalTimingsStr="$finalTimingsStr $finalTiming + "
	    else
		echo "ALERT: Did not find timing information in $finalOutputFile!" >&2
	    fi

	    if [ -n "$nrTimings"  -a  -n "$finalTiming" ]; then
		sumTiming=`echo "$nrTimings $finalTiming" | bc`
	    fi

	    if [ $gzipped -eq 1 ]; then
		rm -f "$finalOutputFile.gz"
		gzip $finalOutputFile
	    fi
	elif [ $DEBUG -eq 1 ]; then
	    echo "Final output file \"$finalOutputFile\" does not exist!" >&2
	fi

	echo "$tag,$dataset$nrIterTiming,$finalTiming,$sumTiming" >> $allTimingOutFile
    done
    if [ "$nrIterTimingsStr" != "3 k 0" ]; then
    	echo "$nrIterTimingsStr p"
	echo "$nrIterTimingsStr p" | dc > $timingOutFile
    elif [ "$finalTimingsStr" != "3 k 0" ]; then
    	echo "$finalTimingsStr p"
	echo "$finalTimingsStr p" | dc > $timingOutFile
    fi
else
    echo "TIMING output file \"$PWD/$allTimingOutFile\" already exists (and I'm not going to recreate it)"
fi
