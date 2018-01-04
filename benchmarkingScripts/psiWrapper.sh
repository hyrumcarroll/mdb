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

# Set MDB_BENCHMARKING_SCRIPTS_DIR to the currently directory unless it is an environmental variable
if [ -z "$MDB_BENCHMARKING_SCRIPTS_DIR" ]; then
    MDB_BENCHMARKING_SCRIPTS_DIR="."
else
    if [ ! -d "$MDB_BENCHMARKING_SCRIPTS_DIR" ]; then
	echo "ERROR: MultiDomainBenchmark benchmarking scripts directory \"$MDB_BENCHMARKING_SCRIPTS_DIR\" does not exist!" >&2
	exit 1
    fi
fi

VERBOSE=1
DEBUG=0
argsStr=""

QUERY=""
FINAL_DATABASE_FULL=""
PSSM_FILE=""
SPOUGE_FILE=""
RELEVANCE_FILE=""
PSI_BLAST_STR="psiblast"
PSI_SEMIGLOBAL_STR="psisemiglobal"
TEST=0
FLAGS_FILE=""
FLAGSFINAL_FILE=""
iterationsDb="$MDB_ITERATIONS_DATABASE"  # default to the environment variable (if set)

removeFiles(){
    while [ "$1" ]; do
	filename="$1"
	shift
	if [ -s "$filename" ]; then
	    cmd="rm -f $filename"
	    if [ $VERBOSE -eq 1 ]; then echo $cmd; fi
	    $cmd
	fi
    done
} # END removeFiles()
    

compressFiles(){
    if [ $# -lt 1 ] ; then
	echo "ERROR: compressFiles() needs arg(s)" >&2
	exit 1
    fi

    while [ "$1" ]; do
	filename="$1"
	shift
	if [ ! -s "$filename" ]; then
	    if [ ! -s "$filename.gz" ]; then
		echo "ERROR: Missing file to compress \"$filename\"!" >&2
		exit 1
	    else
		# a compressed version already exists, next file please
		continue
	    fi
	fi
	if [ -s "$filename.gz" ]; then
	    cmd="rm -f $filename.gz"
	    if [ $DEBUG -eq 1 ]; then echo $cmd; fi
	    $cmd
	fi
	cmd="gzip $filename"
	if [ $DEBUG -eq 1 ]; then echo $cmd; fi
	$cmd
    done
} # END compressFiles()

# gunzip a file specicied by it's decompressed name
# If the gzipped file doesn't exist, but a decompressed does, assume all is good
decompressFiles(){
    if [ $# -lt 1 ] ; then
	echo "ERROR: decompressFiles() needs arg(s)" >&2
	exit 1
    fi

    while [ "$1" ]; do
	filename="$1"
	shift
	if [ ! -s "$filename.gz" ]; then
	    # Compressed files does not exist
	    if [ -s "$filename" ]; then
		# Base does file exists, so don't worry about anything else
		continue
	    else
		echo "ERROR: Missing gzipped version of \"$filename\"!" >&2
		exit 1
	    fi
	fi

	if [ -f "$filename" ]; then
	    cmd="rm -f $filename"
	    if [ $DEBUG -eq 1 ]; then echo $cmd; fi
	    $cmd
	fi
	
	cmd="gunzip $filename.gz"
	if [ $DEBUG -eq 1 ]; then echo $cmd; fi
	$cmd

	# verify something exists
	if [ ! -s $filename ]; then
	    echo "ERROR: Decompressed $filename.gz, but something went wrong...!" >&2
	    exit 1
	fi
    done
} # END decompressFiles()
    
# Determine if the supplied files are complete:
# 1) BLAST output file: Has "Window for multiple hits" (and not "***** No hits found ******", nor "Warning: Failed to retrieve sequence")
# 2) PSSM file: Has "Lambda"
# 3) Blocks output file: Exists (and is not 0-sized)
# Decompress (and leave decompressed) the BLAST output file and the PSSM file
filesComplete(){
    if [ $# -lt 2 ] ; then
	echo "ERROR: filesComplete() needs arg(s)" >&2
	exit 1
    fi

    outputFile="$1"
    pssmCheckpointFile="$2"
    blocksFile="${3-''}"  # if there's a third parameter, use it, otherwise, set it to an empty string
    
    filesComplete=1  # 1: searching complete, 0: searching incomplete

    # Check output file for completeness
    if [ ! -s $outputFile.gz -a ! -s $outputFile ]; then
	return 0
    else
	decompressFiles "$outputFile"
	
	noHitsLine=`grep '***** No hits found ******' $outputFile`
	failedSequenceRetrievalLine=`grep -l 'Warning: Failed to retrieve sequence' $outputFile`
	finishLine=`grep -l 'Window for multiple hits' $outputFile`

	if [ -n "$noHitsLine" ]; then
	    echo "No hits found in $outputFile!" >&2
	    return 0
	elif [ -n "$failedSequenceRetrievalLine" ]; then
	    echo "Redoing iteration search (failed to retrieve sequence in $outputFile)" >&2
	    mv -f $outputFile $outputFile.failed
	    return 0
	elif [ -z "$finishLine" ]; then
	    echo "Redoing iteration search (incomplete output file $outputFile)" >&2
	    mv -f $outputFile $outputFile.failed
	    return 0
	else
	    if [ $VERBOSE -eq 1 ]; then
		convergedLine=`grep -i "Search has CONVERGED!" $outputFile`
		if [ -n "$convergedLine" ]; then
		    echo "Search converged detected in $outputFile";
		fi
	    fi
	    
	    #
	    # Check PSSM file for completeness
	    #
	    if [ ! -s $pssmCheckpointFile.gz -a ! -s $pssmCheckpointFile ]; then
		return 0
	    else
		decompressFiles "$pssmCheckpointFile"
		
		finishLine=`grep -il 'Lambda' $pssmCheckpointFile`
		if [ -z "$finishLine" ]; then
		    echo "Redoing iteration search (incomplete pssmCheckpoint file $pssmCheckpointFile)" >&2
		    mv -f $pssmCheckpointFile $pssmCheckpointFile.failed
		    return 0
		fi
	    fi
	fi	

	if [ "$APP" == "$PSI_SEMIGLOBAL_STR" ]; then
	    # PSI-SemiGLOBAL specific
	    # Test for block file existence
	    if [ ! -s "$blocksFile"  -a  ! -s "$blocksFile.gz" ]; then
		echo "Blocks indices output file \"$blocksFile\" does not exist!"
		return 0
	    fi
	fi
    fi
    
    return 1
    
} # END filesComplete()

usePssmAndScore(){
    #
    # PSI-BLAST/PSI-SemiGLOBAL with saved PSSM on the benchmark database
    #
    
    if [ $# -lt 4 ] ; then
	echo "ERROR: usePssmAndScore() needs arg" >&2
	exit 1
    fi
    
    outputFileLocal="$1"
    hitsOutputFileLocal="$2"
    pssmCheckpointFileLocal="$3"
    SPOUGE_FILELocal="$4"

    runSearch=0
    if [ ! -s $hitsOutputFileLocal.gz -a ! -s $hitsOutputFileLocal ]; then
	runSearch=1
    elif [ ! -s $outputFileLocal.gz -a ! -s $outputFileLocal ]; then
	runSearch=1
    else
	decompressFiles "$hitsOutputFileLocal" "$outputFileLocal"
	
	failedSequenceRetrievalLine=`grep -l 'Warning: Failed to retrieve sequence' $outputFileLocal`
	finishLine=`grep -l 'Window for multiple hits' $hitsOutputFileLocal`
	if [ -n "$failedSequenceRetrievalLine" ]; then
	    if [ $VERBOSE -eq 1 ]; then echo "Redoing final search (Failed to retrieve sequence)"; fi
	    runSearch=1
	    mv -f $outputFileLocal $outputFileLocal.failed
	elif [ -z "$finishLine" ]; then
	    if [ $VERBOSE -eq 1 ]; then echo "Redoing final search"; fi
	    runSearch=1
	    mv -f $hitsOutputFileLocal $hitsOutputFileLocal.failed
	fi
    fi

    if [ $runSearch -eq 1 ]; then
	decompressFiles "$pssmCheckpointFileLocal"
	    
	cmd="$appExecutable -db $finalDatabase  $inputFinalStr  -in_pssm $pssmCheckpointFileLocal  -num_iterations 1 $eValueThresholdFinalStr  -num_threads $NUM_PROCESSORS_TO_USE  -out $hitsOutputFileLocal  $FLAGSFINAL"
	if [ $VERBOSE -eq 1 ]; then echo $cmd; fi
	(echo $cmd; time -p $cmd) &> $outputFileLocal
    
        if [ $? -ne 0 ]; then
    	    echo "ERROR executing \"$cmd\": $!" 1>&2
    	    exit 1;
        fi

	if [ ! -s $hitsOutputFileLocal ]; then
	    echo "ERROR: \"$cmd\" did not produce $hitsOutputFileLocal!" >&2
	    exit 1
	fi
	
        noHitsLine=`grep '***** No hits found ******' $hitsOutputFileLocal`
        if [ -n "$noHitsLine" ]; then
    	    echo "No hits found in $hitsOutputFileLocal! (exiting)" 1>&2
    	    #exit 1;  # don't exit so that any previous iterations can be scored
	    return  
        fi

	failedSequenceRetrievalLine=`grep -l 'Warning: Failed to retrieve sequence' $outputFileLocal`
	if [ -n "$failedSequenceRetrievalLine" ]; then
	    echo "ERROR: Output file has the following error: \"$failedSequenceRetrievalLine\" (final search)!" >&2
	    exit 1
	fi

	finishLine=`grep -l 'Window for multiple hits' $hitsOutputFileLocal`
	if [ -z "$finishLine" ]; then
	    echo "ERROR: hitsOutput file \"$hitsOutputFileLocal\" not complete (final search)!" >&2
	    exit 1
	fi
    else
        if [ $VERBOSE -eq 1 ]; then echo "Skipping final search"; fi
    fi

    
    CLASSIFY_SCRIPT="$MDB_BENCHMARKING_SCRIPTS_DIR/classifyRelevance.pl"
    if [ ! -s $SPOUGE_FILELocal ]; then
	decompressFiles "$hitsOutputFileLocal"

	cmd="$CLASSIFY_SCRIPT   $argsStr  -$APP=$hitsOutputFileLocal  -rel=$RELEVANCE_FILE  --spouge=$SPOUGE_FILELocal --spougeext  $classifyStr"
        if [ $VERBOSE -eq 1 ]; then echo $cmd; fi
        $cmd
    
        if [ $? -ne 0  -o  ! -s $SPOUGE_FILELocal ]; then
	    ls -l $SPOUGE_FILELocal >&2
    	    echo "ERROR executing \"$cmd\": \$?: $?; \$!: $!" 1>&2
    	    exit 1;
        fi
    else
        if [ $VERBOSE -eq 1 ]; then echo "Skipping $CLASSIFY_SCRIPT"; fi
    fi
}  # END usePssmAndScore()


for arg in "$@" ; do

    case "$arg" in
	-verbose     )  VERBOSE=1 ;;
	-v           )  VERBOSE=1 ;;
	-debug       )  DEBUG=1 ;; 
	-d           )  DEBUG=1 ;;

	-psiblast    )  APP="$PSI_BLAST_STR" ;;
	-psisemiglobal   )  APP="$PSI_SEMIGLOBAL_STR" ;;

	-query=*     )  QUERY="`echo $arg | sed -e 's/-query=//'`" ;;
	-msa=*       )  MSA="`echo $arg | sed -e 's/-msa=//'`" ;;

	-iterationsDb=* )  iterationsDb="`echo $arg | sed -e 's/-iterationsDb=//'`" ;;

	-db=*        )  FINAL_DATABASE_FULL="`echo $arg | sed -e 's/-db=//'`" ;;

	-pssm=*      )  PSSM_FILE="`echo $arg | sed -e 's/-pssm=//'`" ;;

	-spouge=*    )  SPOUGE_FILE="`echo $arg | sed -e 's/-spouge=//'`" ;;

	-evalue=*    )  E_VALUE_THRESHOLD_FINAL="`echo $arg | sed -e 's/-evalue=//'`" ;;

	-rel=*       )  RELEVANCE_FILE="`echo $arg | sed -e 's/-rel=//'`" ;;

	-flagsFile=* )  FLAGS_FILE="`echo $arg | sed -e 's/-flagsFile=//'`" ;;
	-flagFile=*  )  FLAGS_FILE="`echo $arg | sed -e 's/-flagFile=//'`" ;;
	-flags=*     )  FLAGS_FILE="`echo $arg | sed -e 's/-flags=//'`" ;;
	-flag=*      )  FLAGS_FILE="`echo $arg | sed -e 's/-flag=//'`" ;;
	-f=*         )  FLAGS_FILE="`echo $arg | sed -e 's/-f=//'`" ;;

	-flagsFinalFile=* )  FLAGSFINAL_FILE="`echo $arg | sed -e 's/-flagsFinalFile=//'`" ;;
	-flagsFinal=*     )  FLAGSFINAL_FILE="`echo $arg | sed -e 's/-flagsFinal=//'`" ;;
		
	-test        )  TEST=1 ;;
	-t           )  TEST=1 ;;

	*            )  if [ -n "$argsStr" ]; then echo "ERROR: this script expects -v and -d AFTER all of args :)!!" !>&2; exit 1; fi ;;
    esac
done 

if [ $DEBUG -eq 1 ]; then
    VERBOSE=1
    argsStr="$argsStr -d"
fi

if [ $VERBOSE -eq 1 ]; then
    argsStr="$argsStr -v"
fi


if [ -n "$FLAGS_FILE"  -a  -s "$FLAGS_FILE" ]; then
    FLAGS=`cat $FLAGS_FILE`
    if [ $DEBUG -eq 1  -o  $VERBOSE -eq 1 ]; then
	echo "DEBUG: FLAGS: $FLAGS"
    fi
fi

if [ -n "$FLAGSFINAL_FILE"  -a  -s "$FLAGSFINAL_FILE" ]; then
    FLAGSFINAL=`cat $FLAGSFINAL_FILE`
    if [ $DEBUG -eq 1  -o  $VERBOSE -eq 1 ]; then
	echo "DEBUG: FLAGSFINAL: $FLAGSFINAL"
    fi
fi

if [ $VERBOSE -eq 1 ]; then 
    echo ""
    echo "$0 args:"
    echo "  APP:                            $APP"
    echo "  QUERY:                          $QUERY"
    echo "  MSA:                            $MSA"
    echo "  iterationsDb:                   $iterationsDb"
    echo "  FINAL_DATABASE_FULL:            $FINAL_DATABASE_FULL"
    echo "  PSSM_FILE:                      $PSSM_FILE"
    echo "  SPOUGE_FILE:                    $SPOUGE_FILE"
    echo "  FLAGS_FILE:                     $FLAGS_FILE"
    echo "  FLAGS:                          $FLAGS"
    echo "  FLAGSFINAL_FILE:                $FLAGSFINAL_FILE"
    echo "  FLAGSFINAL:                     $FLAGSFINAL"
    echo "  TEST:                           $TEST (print commandline, but don't execute runs)"
    echo "  VERBOSE:                        $VERBOSE"
    echo "  DEBUG:                          $DEBUG"
fi

if [ -z "$PSSM_FILE" ]; then
    echo "ERROR: pssm file not specified!" >&2
    exit 1
fi

if [ -s "SPOUGE_FILE" ]; then
    echo "ERROR: spouge file not specified!" >&2
    exit 1
fi


appExecutable="psiblast"  # assume that it's in PATH
if [ "$APP" == "psisemiglobal" ]; then
    appExecutable="$APP"  # assume that it's in PATH
fi

whichTest=`which $appExecutable`
if [ -z "$whichTest" ]; then
    echo "" >&2
    echo "ERROR: Can not find \"$appExecutable\"!" >&2
    echo "" >&2
    exit 1
fi

finalDatabase=${FINAL_DATABASE_FULL%.pin}

if [ -z "$iterationsDb" ]; then
    echo "" >&2
    echo "ERROR: Iterations database not specified (and it must be a full path name)!" >&2
    echo "" >&2
    exit 1
fi    

if [ -n "$E_VALUE_THRESHOLD_FINAL" ]; then
    eValueThresholdFinalStr=" -evalue $E_VALUE_THRESHOLD_FINAL"
else
    eValueThresholdFinalStr=""
fi

#NUM_PROCESSORS_TO_USE=8
NUM_PROCESSORS_TO_USE=1

if [ -z "$QUERY" ]; then
    echo "ERROR: specify -query=<query file>" >&2
    exit 1
fi

if [ ! -s "$QUERY" ]; then
    echo "ERROR: Query file \"$QUERY\" not found: $!" 1>&2
    exit 1;
fi

queryBase=${QUERY##*/} # strip off path
queryBase=${queryBase%.*}

inputIterStr=""
inputFinalStr=""

queryLines=`cat $QUERY | wc -l`
if [ "$queryLines" -gt 2 ]; then
    if [ "$APP" == "$PSI_SEMIGLOBAL_STR" ]; then
	inputIterStr="-in_msa_full $QUERY"
    else
	inputIterStr="-in_msa $QUERY"
    fi
    classifyStr="-family=$queryBase"
else
    inputIterStr="-query $QUERY"
    classifyStr="-taxon=$queryBase"
fi

dbClassifyCriterionFile="$finalDatabase.classifyCriterion"
if [ -s "$dbClassifyCriterionFile" ]; then
    classifyStr="$classifyStr  "`cat $dbClassifyCriterionFile`
else
    classifyStr="$classifyStr  -randomsAsIrrelevants"
fi

echo ""
echo "classifyStr: $classifyStr"
echo ""

#
# PSI-BLAST/PSI-SemiGLOBAL $numTotalIterations iterations on iterationsDb and save checkpoint file (PSSM)
#
pssmCheckpointFile="$PSSM_FILE"
pssmAsciiFile="$pssmCheckpointFile.txt"
outputFile="${queryBase}.$APP$tag.iterationsDb.out"
blocksFile=""  # output from the iteration search, input for the final search
if [ "$APP" == "$PSI_SEMIGLOBAL_STR" ]; then
    blocksFile="${queryBase}${tag}.blocksOut_-iterationsDb.txt"
fi

#
# Look for iterative search files.  If not complet, run iterative search
#
filesComplete $outputFile $pssmCheckpointFile $blocksFile
allFileComplete=$?  # return value from filesComplete

if [ $allFileComplete -eq 0 ]; then
    
    numTotalIterations=5 
    if [ "$APP" == "$PSI_SEMIGLOBAL_STR" ]; then
	# Doesn't the -exit_after_last_pssm_calculation mean that PSI-SemiGLOBAL will stop at the beginning of the last iteration, after the PSSM is calculated?  So, just leave it at the same number of iterations, but the searching on the last is not wasted.
	# # PSI-SemiGLOBAL can save the PSSM and then exit, requiring one less iteration
	# numTotalIterations=$(( numTotalIterations - 1 ))
	
	# Add some flags
	inputIterStr="$inputIterStr  -blocks_out $blocksFile"
	inputIterStr="$inputIterStr  -exit_after_last_pssm_calculation"

	# to avoid any confusion, clear out this output file, if it exists
	removeFiles "$blocksFile.gz"
    fi
    
    cmd="$appExecutable -db $iterationsDb  $inputIterStr  -num_iterations $numTotalIterations -num_threads $NUM_PROCESSORS_TO_USE  -out_pssm $pssmCheckpointFile  -out_ascii_pssm $pssmAsciiFile  $FLAGS"
    if [ $VERBOSE -eq 1 ]; then echo $cmd; fi
    (echo $cmd; time -p $cmd) &> $outputFile

    if [ $? -ne 0 ]; then
	echo "ERROR executing \"$cmd\": \$!: $!, \$?: $?" >&2

	if [ -s $outputFile ]; then
    	    if [ `grep -l "BLAST engine error: No sequences left after purging biased sequences in multiple sequence alignment" $outputFile` ]; then
        	echo "No hits found in $outputFile (after purging biased sequences)" >&2
	    fi
	fi
	exit 1
    fi

    filesComplete $outputFile $pssmCheckpointFile $blocksFile
    allFileComplete=$?  # return value from filesComplete
    if [ $allFileComplete -eq 0 ]; then
	echo "ERROR: Not all files are complete after executing \"$cmd\"" >&2
	exit 1
    fi
else
    if [ $VERBOSE -eq 1 ]; then echo "Skipping iterations search"; fi
fi

filesToCompress="$outputFile $pssmAsciiFile"

#
# Final iteration (searching on a different DB)
#
if [ -z "$SPOUGE_FILE" ]; then
    echo "No -spouge argument, skipping final iteration"
else
    
    outputFileFinal="${queryBase}.$APP$tag.final.out"
    hitsOutputFileFinal="${queryBase}.$APP$tag.hits.final.out"

    if [ "$APP" == "$PSI_SEMIGLOBAL_STR" ]; then
	inputFinalStr=" $inputFinalStr -blocks $blocksFile "
    fi
    
    usePssmAndScore $outputFileFinal  $hitsOutputFileFinal  $pssmCheckpointFile  $SPOUGE_FILE
fi

filesToCompress="$filesToCompress $pssmCheckpointFile $outputFileFinal $hitsOutputFileFinal"
if [ "$APP" == "$PSI_SEMIGLOBAL_STR" ]; then
    filesToCompress="$filesToCompress $blocksFile"
fi

compressFiles $filesToCompress

exit 0;
