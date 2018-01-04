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

#
# If we're doing iterations, then produce a .hmm, otherwise, convert the query into a .hmm.  Then search.
#

# Future: change out --chkhmm for --chkali, and then convert it


# Set MDB_BENCHMARKING_SCRIPTS_DIR to the currently directory unless it is an environmental variable
if [ -z "$MDB_BENCHMARKING_SCRIPTS_DIR" ]; then
    MDB_BENCHMARKING_SCRIPTS_DIR="."
else
    if [ ! -d "$MDB_BENCHMARKING_SCRIPTS_DIR" ]; then
	echo "ERROR: MultiDomainBenchmark benchmarking scripts directory \"$MDB_BENCHMARKING_SCRIPTS_DIR\" does not exist!" >&2
	exit 1
    fi
fi


E_VALUE_THRESHOLD_FINAL=1000

NUM_PROCESSORS_TO_USE=1
ALWAYS_REGENERATE_FILES=0

VERBOSE=1
DEBUG=0

function VERBOSE {
    if [ $VERBOSE -eq 1 ]; then
	echo "$@"
    fi
}

function DEBUG {
    if [ $DEBUG -eq 1 ]; then
	echo "DEBUG: $@"
    fi
}

function run {
    cmd=$1
    shift 1
    products="$@"
    DEBUG "run()"
    DEBUG "    cmd=\"$cmd\""
    DEBUG "    products=\"$products\""
    VERBOSE "$cmd"
    
    eval $cmd
    
    rc=$?
    if [ $rc -ne 0 ]; then
	echo "ERROR: $cmd failed! (rc: $rc)" >&2
	exit 1
    fi
    
    for i in $products; do
	if [ ! -s $i ]; then
	    echo "ERROR: The following script failed to produce (at least) $i: $cmd" >&2
	    exit 1
	fi
    done
}


argsStr=""

#HMMER_NON_ITERATIVE_STR="hmmer-non_iterative"
#HMMER_ITERATIVE_STR="hmmer-iterative"
HMMER_NON_ITERATIVE_STR="phmmer"
HMMER_ITERATIVE_STR="jackhmmer"
QUERY=""
FINAL_DATABASE_FULL_PREFIX=""
SPOUGE_FILE=""
RELEVANCE_FILE=""
TEST=0
FLAGS_FILE=""
FLAGSFINAL_FILE=""
iterationsDb="$MDB_ITERATIONS_DATABASE"  # default to the environment variable (if set)
APP="$HMMER_NON_ITERATIVE_STR" # default
appExecutable="phmmer"         # default


for arg in "$@" ; do
    #echo "arg = $arg"
    case "$arg" in
	-verbose     )  VERBOSE=1 ;;
	-v           )  VERBOSE=1 ;;
	-debug       )  DEBUG=1 ;; 
	-d           )  DEBUG=1 ;;

      # -hmmer-non_iterative|-phmmer  )  APP="$HMMER_NON_ITERATIVE_STR"; appExecutable="$APP" ;;
      # -hmmer-iterative|-jackhmmer   )  APP="$HMMER_ITERATIVE_STR"; appExecutable="$APP" ;;
        -hmmer-non_iterative|-phmmer  )  APP="${arg#-}"; appExecutable="$APP" ;;
        -hmmer-iterative|-jackhmmer   )  APP="${arg#-}"; appExecutable="$APP" ;;

	-query=*     )  QUERY="`echo $arg | sed -e 's/-query=//'`" ;;
	-msa=*       )  MSA="`echo $arg | sed -e 's/-msa=//'`" ;;

	-iterationsDb=* )  iterationsDb="`echo $arg | sed -e 's/-iterationsDb=//'`" ;;

	-db=*        )  FINAL_DATABASE_FULL_PREFIX="`echo $arg | sed -e 's/-db=//'`" ;;

	-spouge=*    )  SPOUGE_FILE="`echo $arg | sed -e 's/-spouge=//'`" ;;

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
    if [ $DEBUG -eq 1 ]; then
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
    echo "  FINAL_DATABASE_FULL_PREFIX:     $FINAL_DATABASE_FULL_PREFIX"
    echo "  SPOUGE_FILE:                    $SPOUGE_FILE"
    echo "  FLAGS_FILE:                     $FLAGS_FILE"
    echo "  FLAGS:                          $FLAGS"
    echo "  FLAGSFINAL_FILE:                $FLAGSFINAL_FILE"
    echo "  FLAGSFINAL:                     $FLAGSFINAL"
    echo "  TEST:                           $TEST (print commandline, but don't execute runs)"
    echo "  VERBOSE:                        $VERBOSE"
    echo "  DEBUG:                          $DEBUG"
fi


if [ -z "$iterationsDb" ]; then
    echo "" >&2
    echo "ERROR: Iterations database not specified (and it must be a full path name)!" >&2
    echo "" >&2
    exit 1
fi
if [ "$iterationsDb" == "${iterationsDb%.fa}" ]; then
    echo "NOTE: Appending .fa to the end of Iterations DB filename \"$iterationsDb\"" >&2
    iterationsDb="$iterationsDb.fa"
fi

for ext in pin fa; do
    if [ "$FINAL_DATABASE_FULL_PREFIX" != "${FINAL_DATABASE_FULL_PREFIX%.$ext}" ]; then
	# strip off extension
	FINAL_DATABASE_FULL_PREFIX="${FINAL_DATABASE_FULL_PREFIX%.$ext}"
    fi
done

finalDatabase="$FINAL_DATABASE_FULL_PREFIX.fa"
if [ ! -s "$finalDatabase" ]; then
    echo "ERROR: Final database \"$finalDatabase\" does not exist!" >&2
    exit 1
fi

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

# inputStr=""
# singleSeqQuery=1
queryLines=`cat $QUERY | wc -l`
if [ "$queryLines" -gt 2 ]; then
    # singleSeqQuery=0
    # inputStr="-in_msa $QUERY"
    classifyStr="-family=$queryBase"
else
    # singleSeqQuery=1
    # inputStr="-query $QUERY"
    classifyStr="-taxon=$queryBase"
fi

dbClassifyCriterionFile="$FINAL_DATABASE_FULL_PREFIX.classifyCriterion"
if [ -s "$dbClassifyCriterionFile" ]; then
    if [ $DEBUG -eq 1 ]; then echo "Reading dbClassifyCriterionFile \"$dbClassifyCriterionFile\" . . ." >&2; fi
    classifyStr="$classifyStr  "`cat $dbClassifyCriterionFile`
else
    classifyStr="$classifyStr  -randomsAsIrrelevants"
fi

# append command-line classification args
classifyStr="$classifyStr  $CLASSIFY_ARGS"

if [ $VERBOSE -eq 1 ]; then
    echo "classifyStr: $classifyStr"
fi

hmmbuild="hmmbuild"
whichTest=`which $hmmbuild`
if [ -z "$whichTest" ]; then
    echo "" >&2
    echo "ERROR: Can not find \"$hmmbuild\"!" >&2
    echo "" >&2
fi

hmmsearch="hmmsearch"
whichTest=`which $hmmsearch`
if [ -z "$whichTest" ]; then
    echo "" >&2
    echo "ERROR: Can not find \"$hmmsearch\"!" >&2
    echo "" >&2
fi


numTotalIterations=0
if [ "$APP" == "$HMMER_ITERATIVE_STR" ]; then
    numTotalIterations=5
else
    numTotalIterations=0
fi

tag=""
finalHmmInputFilename="${queryBase}.$APP$tag.finalIn.hmm"
generateFinalHmmInput=0
if [ ! -s "$finalHmmInputFilename" ]; then
    generateFinalHmmInput=1
else
    finishLine=`tail -n 1 $finalHmmInputFilename | grep '^//$'`
    if [ -z "$finishLine" ]; then
	generateFinalHmmInput=1
    fi
fi
if [ $generateFinalHmmInput -eq 1 ]; then 
    if [ $numTotalIterations -eq 0 ]; then
	#
	# Non-iterative
	#
	run "$hmmbuild $finalHmmInputFilename $QUERY"  "$finalHmmInputFilename"
    else # iterations -eq 0
	#
	# Iterative searching
	#
	
	tag=""

	#
	# $numTotalIterations iterations on iterationsDb database and save checkpoint MSA
	#
	outputFile="${queryBase}.$APP$tag.iterationsDb.out"
	hmmCheckpointFilePrefix="${queryBase}.$APP$tag.iterationsDb" # At the start of each iteration, checkpoint the query HMM, saving it to a file named <prefix>-<n>.hmm where <n> is the iteration number (from 1..N).

	runSearch=$ALWAYS_REGENERATE_FILES
	if [ ! -s $finalHmmInputFilename.gz -a ! -s $finalHmmInputFilename  -o  ! -s $outputFile.gz -a ! -s $outputFile ]; then
	    runSearch=1
	else
	    if [ -s $finalHmmInputFilename.gz -o -s $finalHmmInputFilename ]; then
		if [ ! -s $finalHmmInputFilename -a -s $finalHmmInputFilename.gz ]; then
		    run "rm -f $finalHmmInputFilename;  gunzip $finalHmmInputFilename.gz"
		fi
		finishLine=`tail -n 1 $finalHmmInputFilename | grep '^//$'`
		if [ -z "$finishLine" ]; then
		    if [ $VERBOSE -eq 1 ]; then echo "Redoing initial search (numTotalIterations: $numTotalIterations) ($finalHmmInputFilename is incomplete)"; fi
		    runSearch=1
		    run "rm \"$finalHmmInputFilename\""
		    # else
		    # leave the .hmm file unzipped
		    # 	rm -f $finalHmmInputFilename.gz;   gzip $finalHmmInputFilename	    
		fi
	    fi
	    if [ -s $outputFile.gz -o -s $outputFile ]; then
		if [ ! -s $outputFile -a -s $outputFile.gz ]; then
		    run "rm -f $outputFile;  gunzip $outputFile.gz"
		fi
		finishLine=`grep -l '^\[ok\]$' $outputFile`
		if [ -z "$finishLine" ]; then
		    if [ $VERBOSE -eq 1 ]; then echo "Redoing initial search (numTotalIterations: $numTotalIterations) ($outputFile is incomplete)"; fi
		    runSearch=1
		else
		    run "rm -f $outputFile.gz;   gzip $outputFile"
		fi
	    fi
	fi

	if [ $runSearch -eq 1 ]; then
	    # From the HMMer Userguide:
	    # --noali Omit the alignment section from the main output. This can greatly reduce the output volume. [Alignments in the output are necessary to calculate the alignment overlap]
	    # -N <n> Set the maximum number of iterations to <n>. The default is 5. If N=1, the result is equivalent to a phmmer search.
	    # -A <f> After the final iteration, save an annotated multiple alignment of all hits satisfying inclusion thresholds (also including the original query) to <f> in Stockholm format.
	    # --chkali <prefix> At the end of each iteration, checkpoint an alignment of all domains satisfying inclusion thresholds (e.g. what will become the query HMM for the next iteration), saving it to a file named <checkpoint file prefix>-<n>.sto in Stockholm format, where <n> is the iteration number (from 1..N).
	    # --chkhmm <prefix> Direct the human-readable output to a file <f>. At the start of each iteration, checkpoint the query HMM, saving it to a file named <prefix>-<n>.hmm where <n> is the iteration number (from 1..N).
	    
    	    run "(time -p $appExecutable  $FLAGS --noali -N $numTotalIterations  --cpu $NUM_PROCESSORS_TO_USE  --chkhmm $hmmCheckpointFilePrefix  $QUERY  $iterationsDb | grep '^\[ok\]$') &> $outputFile"  # "$finalIterationInputFilename" #  -A $finalIterationInputFilename --chkali $alignmentCheckpointFilenamePrefilx # --incE $E_VALUE_THRESHOLD_INITIAL 
	    #
            # Need to update for HMMer 	
	    #
	    # if [ -n "$noHitsLine" ]; then
	    #     echo "No hits found in $outputFile! (skipping)" 1>&2
	    #     continue
	    # fi

	    finishLine=`grep -l '^\[ok\]$' $outputFile`
	    if [ -z "$finishLine" ]; then
		echo "ERROR: output file \"$outputFile\" not complete (initial search)!" >&2
		exit 1
	    fi

	    # make a link from the last .hmm iteration file to the final iteration input file
	    lastSuccessfulIter=$numTotalIterations
	    while [ $lastSuccessfulIter -gt 0 ]; do
		if [ -s "$hmmCheckpointFilePrefix-$lastSuccessfulIter.hmm" ]; then
		    break
		elif [  -s "$hmmCheckpointFilePrefix-$lastSuccessfulIter.hmm.gz" ]; then
		    run "rm -f \"$hmmCheckpointFilePrefix-$lastSuccessfulIter.hmm\";  gunzip \"$hmmCheckpointFilePrefix-$lastSuccessfulIter.hmm.gz\""
		    break
		fi
		lastSuccessfulIter=$(( lastSuccessfulIter - 1 ))
	    done
	    
	    lastIterationHmmFilename="$hmmCheckpointFilePrefix-$lastSuccessfulIter.hmm"
	    if [ ! -s "$lastIterationHmmFilename" ]; then
		echo "ERROR: Last iteration .hmm file \"$lastIterationHmmFilename\" not found!" >&2
		exit 1
	    else
		run "ln -s \"$lastIterationHmmFilename\" \"$finalHmmInputFilename\""
	    fi
	    
	    run "rm -f $outputFile.gz;   gzip $outputFile"
	    
	else
	    if [ $VERBOSE -eq 1 ]; then echo "Skipping initial search (numTotalIterations: $numTotalIterations)"; fi
	fi

	if [ ! -s $outputFile.gz ]; then
	    if [ ! -s $outputFile ]; then
		ls -l $outputFile $outputFile.gz >&2
		echo "ERROR: $outputFile.gz: $! (numTotalIterations: $numTotalIterations)" >&2
		exit 1
	    else
		run "rm -f $outputFile.gz;   gzip $outputFile"
	    fi
	fi

	# for i in $hmmCheckpointFilePrefix-*; do
	#     if [ -s "$i" ]; then
	# 	gzip "$i"
	#     fi
	# done
    fi  # iterations -eq 0

    # validate finalHmmInputFilename
    if [ ! -s "$finalHmmInputFilename" ]; then
	echo "ERROR: HMM input file \"$finalHmmInputFilename\" does not exist!" >&2
	exit 1
    else
	finishLine=`tail -n 1 $finalHmmInputFilename | grep '^//$'`
	if [ -z "$finishLine" ]; then
	    echo "ERROR: HMM input file \"$finalHmmInputFilename\" not complete!" >&2
	    exit 1
	fi
    fi
fi

#
# Search the benchmark database with 1) HMM generated from a single sequence query or 2) the last iteration HMM 
#
hitsOutputFile="${queryBase}.$APP$tag.hits.final.out"
outputFile="${queryBase}.$APP$tag.final.out"
runSearch=$ALWAYS_REGENERATE_FILES
if [ ! -s $hitsOutputFile.gz -a ! -s $hitsOutputFile ]; then
    runSearch=1
# elif [ ! -s $outputFile.gz -a ! -s $outputFile ]; then
#     runSearch=1
else
    if [ ! -s $hitsOutputFile -a -s $hitsOutputFile.gz ]; then
	run "rm -f $hitsOutputFile;  gunzip $hitsOutputFile.gz"
    fi
    finishLine=`grep -l '^\[ok\]$' $hitsOutputFile`
    if [ -z "$finishLine" ]; then
	if [ $VERBOSE -eq 1 ]; then echo "Redoing final search (numTotalIterations: $numTotalIterations)"; fi
	runSearch=1
    else
	run "rm -f $hitsOutputFile.gz;   gzip $hitsOutputFile"
    fi
fi

if [ $runSearch -eq 1 ]; then
    if [ ! -s $finalHmmInputFilename ]; then
	if [ ! -s $finalHmmInputFilename.gz ]; then
	    echo "ERROR: $finalHmmInputFilename: $!" >&2
	    exit 1;
	else
	    run "rm -f $finalHmmInputFilename;   gunzip $finalHmmInputFilename.gz"
	fi
    fi

    # run "time -p hmmsearch $FLAGS -E $E_VALUE_THRESHOLD_FINAL --cpu $NUM_PROCESSORS_TO_USE -o $hitsOutputFile  $finalHmmInputFilename  $finalDatabase" # &> $outputFile" #  --tblout $hitsOutputFile  --noali
    run "(time -p hmmsearch $FLAGS -E $E_VALUE_THRESHOLD_FINAL --cpu $NUM_PROCESSORS_TO_USE -o $hitsOutputFile  $finalHmmInputFilename  $finalDatabase) &> $outputFile" #  --tblout $hitsOutputFile  --noali

    finishLine=`grep -l '^\[ok\]$' $hitsOutputFile`
    if [ -z "$finishLine" ]; then
	echo "ERROR: Hits output file \"$hitsOutputFile\" not complete (final search)!" >&2
	exit 1
    fi

    # if [ $numTotalIterations -ne 0 ]; then
    # 	run "rm -f $finalHmmInputFilename.gz;   gzip $finalHmmInputFilename"
    # fi
else
    if [ $VERBOSE -eq 1 ]; then echo "Skipping final search (numTotalIterations: $numTotalIterations)"; fi
fi


CLASSIFY_SCRIPT="$MDB_BENCHMARKING_SCRIPTS_DIR/classifyRelevance.pl"
if [  ! -s $SPOUGE_FILE ]; then

    if [ ! -s $hitsOutputFile ]; then
	run "rm -f $hitsOutputFile;   gunzip $hitsOutputFile.gz"
    fi

    cmd="$CLASSIFY_SCRIPT   $argsStr  -$APP=$hitsOutputFile  -rel=$RELEVANCE_FILE  --spouge=$SPOUGE_FILE  --spougeext $classifyStr"
    if [ $VERBOSE -eq 1 ]; then echo $cmd; fi
    $cmd
    
    if [ $? -ne 0  -o  ! -s $SPOUGE_FILE ]; then
	ls -l $SPOUGE_FILE >&2
    	echo "ERROR executing \"$cmd\": \$?: $?; \$!: $!" 1>&2
    	exit 1;
    fi

    run "rm -f $hitsOutputFile.gz;   gzip $hitsOutputFile"
else
    if [ $VERBOSE -eq 1 ]; then echo "Skipping $CLASSIFY_SCRIPT  (numTotalIterations: $numTotalIterations)"; fi
fi

if [ ! -s $hitsOutputFile.gz ]; then
    if [ ! -s $hitsOutputFile ]; then
	ls -l $hitsOutputFile $hitsOutputFile.gz >&2
	echo "ERROR: $hitsOutputFile.gz: $! (numTotalIterations: $numTotalIterations)" >&2
	exit 1
    else
	run "rm -f $hitsOutputFile.gz;   gzip $hitsOutputFile"
    fi
fi

exit 0;
