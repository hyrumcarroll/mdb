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

DEBUG=0
#DEBUG=1

# If MDB_DB_BASE_DIR is not an environmental variable, then set it to the currect directory
if [ -z "$MDB_DB_BASE_DIR" ]; then
    MDB_DB_BASE_DIR="."
fi

#
# subroutines
#

exitIfFileDoesNotExist(){
    if [ ! -s $1 ]; then
	if [ -n "$2" ]; then 
	    errorStr="$2"
	else
	    errorStr="ERROR: \"$1\" does not exist (or is zero length)"
	fi
	echo "$errorStr" >&2
	    
	exit 1
    fi
}


# process command-line args


if [ -z "$1" ]; then
    echo "Usage: $0 <directory of directories that have retrieval results> <TAP-k output file base name> <relevance file>" >&2
    exit 1;
fi

runsDir=$1
shift
if [ $DEBUG -eq 1 ]; then echo "DEBUG: runsDir: $runsDir"; fi

tapOutBase="$1"
shift
if [ $DEBUG -eq 1 ]; then echo "DEBUG: tapOutBase: $tapOutBase"; fi

relevanceFile=$1
shift
if [ $DEBUG -eq 1 ]; then echo "DEBUG: relevanceFile: $relevanceFile"; fi

# performance optimization
eValueThreshold=$1
shift
if [ "$eValueThreshold" = "." ]; then
    eValueThreshold=""
    thresholdStr=""
else
    thresholdStr=".$eValueThreshold"
fi
if [ $DEBUG -eq 1 ]; then echo "DEBUG: eValueThreshold: $eValueThreshold"; fi
thresholdStr=".$eValueThreshold"

#queryType
queryType=$1
shift
if [ $DEBUG -eq 1 ]; then echo "DEBUG: queryType: $queryType"; fi

benchmarkDB=$1
shift
if [ $DEBUG -eq 1 ]; then echo "DEBUG: benchmarkDB: $benchmarkDB"; fi
if [ ! -d "$MDB_DB_BASE_DIR/$benchmarkDB/" ]; then
    echo "ERROR: Directory \"$MDB_DB_BASE_DIR/$benchmarkDB/\" (with specified benchmarkDB: $benchmarkDB) does not exist!" >&2
    exit 1
fi

queriesFile="$MDB_DB_BASE_DIR/$benchmarkDB/queriesList.txt"
exitIfFileDoesNotExist $queriesFile
QUERIES=`cat $queriesFile`
if [ $DEBUG -gt 1 ]; then echo "QUERIES: $QUERIES"; fi



spougeFileDir="${runsDir%/}"
spougeFileDir="${spougeFileDir%/*}"

cmd="cd $spougeFileDir"
echo $cmd
$cmd

ks="1 3 5 20"
tapOutFile="$tapOutBase.tapKs.out"
if [ ! -s "$tapOutFile" ]; then
    submitDir=`basename $runsDir`
    spougeFile="${submitDir}-all.spouge"

    if [ $DEBUG -eq 1 ]; then echo "DEBUG: spougeFile: $spougeFile (spougeFileDir: $spougeFileDir)"; fi

    if [ ! -s "$spougeFile" ]; then
	if [ $DEBUG -eq 1 ]; then echo "DEBUG: Listing .spouge files . . ."; fi
	retrevialFileSpouges=""

	for queryFile in $QUERIES; do
	    if [ $DEBUG -gt 1 ]; then echo "queryFile: $queryFile"; fi
            queryFileBase=${queryFile##*/}
            queryFileBase=${queryFileBase%.*}

	    dir="$submitDir/$queryFileBase.$tapOutBase"
	#
	# for dir in $submitDir/*; do
	#
	    if [ ! -d "$dir" ]; then
		echo "ERROR: Expecting \"$dir\" to be a directory!" >&2
		exit 1
	    fi

	    baseDir=`basename $dir`
    	    retrevialFileSpougeBase="$dir/$baseDir.spouge"
    	    retrevialFileSpouge="$retrevialFileSpougeBase$thresholdStr"
    	    if [ ! -s "$retrevialFileSpouge" ]; then
    		if [ -n "$eValueThreshold"  -a   -s "$retrevialFileSpougeBase" ]; then
		    if [ ! -s "$retrevialFileSpouge" ]; then
		        # truncate .spouge file
    			cmd="./spouge2spougeE.pl -spouge=$retrevialFileSpougeBase -spougee=$retrevialFileSpouge -spougeeEValue=$eValueThreshold"
			if [ $DEBUG -eq 1 ]; then echo $cmd; fi
			$cmd
		    fi
		else
		    retrevialFile="$dir/$baseDir.tab"
    		    if [ -s "$retrevialFile" ]; then
			echo "ERROR: Should not get here!" >&2
			exit 1
		    else
			hitsFile="hits_$baseDir.final.out"
			if [ -s "$hitsFile"  -o  -s "$hitsFile.gz" ]; then
			    echo "Retrieval file, $retrevialFile, not found, but the hits file was ($hitsFile)!" >&2
			    exit 1
			fi 
			family="${baseDir%%\.*}"
			echo "$retrevialFile" > $retrevialFileSpouge
			cmd="./numRelevantRecords.pl -rel=$relevanceFile -${queryType}=$family  -ignoreSelfHit"
			echo $cmd
			$cmd >> $retrevialFileSpouge
			echo "" >> $retrevialFileSpouge
		    fi
    		fi
	    fi
	    retrevialFileSpouges="$retrevialFileSpouges $retrevialFileSpouge"
	done
	
	for i in $retrevialFileSpouges; do cat $i; done > $spougeFile

	# hack
	# for some reason, cat in the loop above is truncating trailing newlines
	perl -pi -e 's/\n+/~/g' $spougeFile
	echo "" >> $spougeFile
	perl -pi -e 's/~(\D)/\n\n$1/g' $spougeFile
	perl -pi -e 's/~/\n/g' $spougeFile
    fi

    tmpOutFile=".$$.out"
    cmd="tap.pl -i $spougeFile -k $ks -w $tapOutFile"
    echo $cmd
    $cmd &> $tmpOutFile
    rc=$?
    cat $tmpOutFile
    
    if [ $rc -ne 0 ]; then
	# Fewer than 0.5 of the retrieval lists have 3 errors.
    	if [ `grep -l "Fewer than " $tmpOutFile` ]; then
	    errorK=`grep "Fewer than" $tmpOutFile`
	    errorK=${errorK% errors.}
	    errorK=${errorK##* }
	    echo "ALERT: Calculating TAP-ks up to $errorK"
	    tapOutFiles=""
	    for k in $ks; do
		if [ $k -lt $errorK ]; then
		    cmd="tap.pl -i $spougeFile -k $k -w $tapOutFile.$k"
		    echo $cmd
		    $cmd
		    tapOutFiles="$tapOutFiles $tapOutFile.$k"
		fi
	    done
	    if [ -n "$tapOutFiles" ]; then 
		cat $tapOutFiles > $tapOutFile
	    fi
	fi
    fi
    rm $tmpOutFile
else
    echo "TAP output file already exists (and I'm not going to recreate it)"
fi

for tapK in $ks; do
    tapKFile="$tapOutBase.tap${tapK}"
    grep -A 1 "^EPQ" $tapOutFile | grep "^$tapK" > $tapKFile
    perl -pi -e 's/^\d+\s+\([\d\.e\-]+\)\s+([\d\.e\-]+)/$1/' $tapKFile
    if [ $DEBUG -eq 1 ]; then echo "DEBUG: TAP-$tapK value: "`cat $tapKFile`; fi
done
