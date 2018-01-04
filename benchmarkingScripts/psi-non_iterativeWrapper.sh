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

ALWAYS_REGENERATE_FILES=0

VERBOSE=1
DEBUG=0
argsStr=""

QUERY=""
FINAL_DATABASE_FULL=""
TAB_FILE=""
SPOUGE_FILE=""
E_VALUE_THRESHOLD=""
RELEVANCE_FILE=""
PSI_BLAST_NON_ITERATIVE_STR="psiblast-non_iterative"
PSI_BLAST_STR="psiblast"
PSI_SEMIGLOBAL_NON_ITERATIVE_STR="psisemiglobal-non_iterative"
#BLASTP_STR="blastp"
PSI_SEMIGLOBAL_STR="psisemiglobal"
TEST=0
FLAGS_FILE=""
APP="$PSI_BLAST_NON_ITERATIVE_STR" # defaults to psiblast
CLASSIFY_HOOK=""
CLASSIFY_ARGS=""

for arg in "$@" ; do
    #echo "arg = $arg"
    case "$arg" in
	-verbose     )  VERBOSE=1 ;;
	-v           )  VERBOSE=1 ;;
	-debug       )  DEBUG=1 ;; 
	-d           )  DEBUG=1 ;;

	-psiblast-non_iterative ) APP="$PSI_BLAST_NON_ITERATIVE_STR" ;;
	-psisemiglobal-non_iterative   )  APP="$PSI_SEMIGLOBAL_NON_ITERATIVE_STR" ;;
	-psisemiglobal   )  APP="$PSI_SEMIGLOBAL_STR" ;;

	-query=*     )  QUERY="`echo $arg | sed -e 's/-query=//'`" ;;
	-msa=*       )  MSA="`echo $arg | sed -e 's/-msa=//'`" ;;

	-db=*        )  FINAL_DATABASE_FULL="`echo $arg | sed -e 's/-db=//'`" ;;

	-out=*       )  OUT_FILE="`echo $arg | sed -e 's/-out=//'`" ;;
	
	-hits=*      )  HITS_FILE="`echo $arg | sed -e 's/-hits=//'`" ;;
	
	-tab=*       )  TAB_FILE="`echo $arg | sed -e 's/-tab=//'`" ;;

	-spouge=*    )  SPOUGE_FILE="`echo $arg | sed -e 's/-spouge=//'`" ;;

	-evalue=*    )  E_VALUE_THRESHOLD="`echo $arg | sed -e 's/-evalue=//'`" ;;

	-roc50=*     )  ROC_FILE="`echo $arg | sed -e 's/-spouge=//'`" ;;

	-rel=*       )  RELEVANCE_FILE="`echo $arg | sed -e 's/-rel=//'`" ;;

	-flagsFile=* )  FLAGS_FILE="`echo $arg | sed -e 's/-flagsFile=//'`" ;;
	-flagFile=*  )  FLAGS_FILE="`echo $arg | sed -e 's/-flagFile=//'`" ;;
	-flags=*     )  FLAGS_FILE="`echo $arg | sed -e 's/-flags=//'`" ;;
	-flag=*      )  FLAGS_FILE="`echo $arg | sed -e 's/-flag=//'`" ;;
	-f=*         )  FLAGS_FILE="`echo $arg | sed -e 's/-f=//'`" ;;

	-classifyHook=*)  CLASSIFY_HOOK="`echo $arg | sed -e 's/-classifyHook=//'`" ;;
	-classifyArgs=*)  CLASSIFY_ARGS=" $CLASSIFY_ARGS `echo $arg | sed -e 's/-classifyArgs=//'`" ;;
	
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

if [ $VERBOSE -eq 1 ]; then 
    echo ""
    echo "$0 args:"
    echo "  APP:                            $APP"
    echo "  QUERY:                          $QUERY"
    echo "  MSA:                            $MSA"
    echo "  FINAL_DATABASE_FULL:            $FINAL_DATABASE_FULL"
    echo "  TAB_FILE:                       $TAB_FILE"
    echo "  SPOUGE_FILE:                    $SPOUGE_FILE"
    echo "  E_VALUE_THRESHOLD:              $E_VALUE_THRESHOLD"
    echo "  FLAGS_FILE:                     $FLAGS_FILE"
    echo "  FLAGS:                          $FLAGS"
    echo "  CLASSIFY_HOOK:                  $CLASSIFY_HOOK"
    echo "  CLASSIFY_ARGS:                  $CLASSIFY_ARGS"
    echo "  TEST:                           $TEST (print commandline, but don't execute runs) (NYI!)"
    echo "  VERBOSE:                        $VERBOSE"
    echo "  DEBUG:                          $DEBUG"
fi

baseApp="${APP%-non_iterative}" # strip off "-non_iterative" if present

appExecutable="psiblast"  # assume that it's in PATH
if [ "$APP" == "$PSI_SEMIGLOBAL_STR"  -o  "$APP" == "$PSI_SEMIGLOBAL_NON_ITERATIVE_STR" ]; then
    appExecutable="psisemiglobal"  # assume that it's in PATH
fi

whichTest=`which $appExecutable`
if [ -z "$whichTest" ]; then
    echo "" >&2
    echo "ERROR: Can not find \"$appExecutable\"!" >&2
    echo "" >&2
    exit 1
fi

finalDatabase=${FINAL_DATABASE_FULL%.pin}

if [ -n "$E_VALUE_THRESHOLD" ]; then
    eValueThresholdFinalStr=" -evalue $E_VALUE_THRESHOLD"
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

inputStr=""
inputFinalStr=""

queryLines=`cat $QUERY | wc -l`
if [ $DEBUG -eq 1 ]; then echo "queryLines: $queryLines" >&2; fi
if [ "$queryLines" -gt 2 ]; then
    if [ "$baseApp" == "$PSI_SEMIGLOBAL_STR" ]; then
	inputStr="-in_msa_full $QUERY"
    else
	inputStr="-in_msa $QUERY"
    fi
    classifyStr="-family=$queryBase"
else
    inputStr="-query $QUERY"
    classifyStr="-taxon=$queryBase"
fi

dbClassifyCriterionFile="$finalDatabase.classifyCriterion"
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

blocksOutputFile="blocksOut_${queryBase}$tag-nr.txt"
if [ "$baseApp" == "$PSI_SEMIGLOBAL_STR" ]; then
    inputStr="$inputStr -blocks_out $blocksOutputFile"
fi

tag=""

if [ -n "$OUT_FILE" ]; then
    outputFile="$OUT_FILE"
else
    outputFile="${queryBase}.$APP$tag.final.out"
fi
if [ -n "$HITS_FILE" ]; then
    hitsOutputFile="$HITS_FILE"
else
    hitsOutputFile="${queryBase}.$APP$tag.hits.final.out"
fi

if [ ! -s $hitsOutputFile.gz -a ! -s $hitsOutputFile ]; then
    runSearch=1
else
    runSearch=$ALWAYS_REGENERATE_FILES
    if [ ! -s $hitsOutputFile -a -s $hitsOutputFile.gz ]; then
	rm -f $hitsOutputFile;  gunzip $hitsOutputFile.gz
    fi
    finishLine=`grep -l 'Window for multiple hits' $hitsOutputFile`
    if [ -z "$finishLine" ]; then
	if [ $VERBOSE -eq 1 ]; then echo "Redoing final search"; fi
	runSearch=1
    else
	rm -f $hitsOutputFile.gz;   gzip $hitsOutputFile	    
    fi
fi

if [ $runSearch -eq 1 ]; then
    cmd="$appExecutable -db $finalDatabase  $inputStr  -num_threads $NUM_PROCESSORS_TO_USE   $eValueThresholdFinalStr   -out $hitsOutputFile  $FLAGS "
    if [ $VERBOSE -eq 1 ]; then echo $cmd; fi
    (echo $cmd; time -p $cmd) &> $outputFile
    
    if [ $? -ne 0 ]; then
    	echo "ERROR executing \"$cmd\": $!" 1>&2
    	exit 1;
    fi
    
    noHitsLine=`grep '***** No hits found ******' $hitsOutputFile`
    if [ -n "$noHitsLine" ]; then
    	echo "No hits found in $hitsOutputFile! (exiting)" 1>&2
    	exit 1;
    fi

    finishLine=`grep -l 'Window for multiple hits' $hitsOutputFile`
    if [ -z "$finishLine" ]; then
	echo "ERROR: hitsOutput file \"$hitsOutputFile\" not complete (final search)!" >&2
	exit 1
    fi

    rm -f $outputFile.gz;   gzip $outputFile
else
    if [ $VERBOSE -eq 1 ]; then echo "Skipping final search"; fi
fi


CLASSIFY_SCRIPT="$MDB_BENCHMARKING_SCRIPTS_DIR/classifyRelevance.pl"

classify=$ALWAYS_REGENERATE_FILES
if [ -n "$SPOUGE_FILE"  -a  ! -s "$SPOUGE_FILE" ]; then
    classify=1
fi

if [ $DEBUG -eq 1 ]; then echo "classify: $classify (SPOUGE_FILE: $SPOUGE_FILE)" >&2; fi

if [ $classify -eq 1 ]; then

    if [ ! -s $hitsOutputFile ]; then
	rm -f $hitsOutputFile;   gunzip $hitsOutputFile.gz
    fi
    

    if [ -n "$CLASSIFY_HOOK" ]; then
	if [ $DEBUG -eq 1 ]; then echo "classify hook: $CLASSIFY_HOOK" >&2; fi
	$CLASSIFY_HOOK
    fi
    

    cmd="$CLASSIFY_SCRIPT   $argsStr  -$baseApp=$hitsOutputFile  -rel=$RELEVANCE_FILE  --spouge=$SPOUGE_FILE  --spougeext  $classifyStr"
    if [ $VERBOSE -eq 1 ]; then echo $cmd; fi
    $cmd
    
    if [ $? -ne 0  -o  ! -s $SPOUGE_FILE ]; then
	ls -l $SPOUGE_FILE >&2
    	echo "ERROR executing \"$cmd\": \$?: $?; \$!: $!" 1>&2
    	exit 1;
    fi

    rm -f $hitsOutputFile.gz;   gzip $hitsOutputFile
else
    if [ $VERBOSE -eq 1 ]; then echo "Skipping $CLASSIFY_SCRIPT "; fi
fi

if [ ! -s $hitsOutputFile.gz ]; then
    if [ ! -s $hitsOutputFile ]; then
	ls -l $hitsOutputFile $hitsOutputFile.gz >&2
	echo "ERROR: $hitsOutputFile.gz: $!" >&2
	exit 1
    else
	rm -f $hitsOutputFile.gz;   gzip $hitsOutputFile
    fi
fi
    
exit 0;
