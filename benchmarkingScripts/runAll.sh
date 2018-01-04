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

DO_NON_ITERATIVE=0
DO_ITERATIVE=0

MULTIDOMAINBENCHMARK_TRAINING=0
MULTIDOMAINBENCHMARK_TEST=0

#pssmOnlyStr="-pssmOnly"
pssmOnlyStr=""

VERBOSE=1
DEBUG=0

BLAST_COMMON_ARGS=" -num_descriptions 9999  -num_alignments 9999 "  # necessary for ROC_n analysis
EVALUE_STR="-evalue=1000"


usageStr="Usage $0:
	-noniterative | -non-iterative | -non_iterative 
	-iterative
	
	-training | -MultiDomainBenchmark-training
	-test |     -MultiDomainBenchmark-test
	
	[-queuecount=* | -queue=*]
	[-totalcount=* | -total=*]
	
	[-v | -verbose]
	[-nov | -noverbose]
	[-d | -debug]
		          
	[-qsub]
	[-sh]
	[-resubmit]
	[-t | -testOnly]
"         

for arg in "$@" ; do
    #echo "arg = $arg"
    case "$arg" in
	-v|-verbose       )  VERBOSE=1 ;;
	-nov|-noverbose   )  VERBOSE=0 ;;
	-d|-debug         )  DEBUG=1 ;; 
		          
	-qsub             )  executeStr="$arg" ;;
	-sh               )  executeStr="$arg" ;;
	-resubmit         )  executeStr="-qsub $arg" ;;
	-t | -testOnly | -check   )  executeStr="" ;;
		          
	# -count=*          )  countStr="$arg" ;;
        # -c=*              )  countStr="-count="`echo $arg | sed -e 's/-c=//'` ;;		
	-queuecount=*     )  countStr="$arg" ;;
        -queue=*          )  countStr="$arg" ;;
	-totalcount=*     )  countStr="$arg" ;;
        -total=*          )  countStr="$arg" ;;
		          
	-noniterative | -non-iterative | -non_iterative  )  DO_NON_ITERATIVE=1 ;;
	-iterative        )  DO_ITERATIVE=1 ;;
		          
	# DBs	          
	-training | -MultiDomainBenchmark-training  )  MULTIDOMAINBENCHMARK_TRAINING=1 ;;
	-test |     -MultiDomainBenchmark-test      )  MULTIDOMAINBENCHMARK_TEST=1 ;;
		       
	*                 )  echo "ERROR: Unknown arg: $arg!" >&2; exit 1;;
    esac
done 


verbosityStr=""
if [ $VERBOSE -gt 0 ]; then
    verbosityStr="$verbosityStr -v"
fi
if [ $DEBUG -gt 0 ]; then
    verbosityStr="$verbosityStr -d"
fi

if [ $DEBUG -eq 1 ]; then 
    echo "DO_NON_ITERATIVE:   $DO_NON_ITERATIVE"
    echo "DO_ITERATIVE:       $DO_ITERATIVE"
    echo "MULTIDOMAINBENCHMARK_TRAINING:  $MULTIDOMAINBENCHMARK_TRAINING"
    echo "MULTIDOMAINBENCHMARK_TEST:      $MULTIDOMAINBENCHMARK_TEST"
    echo "executeStr: $executeStr"
    echo "countStr: $countStr"
    echo "verbosityStr: $verbosityStr"
    echo ""
fi

if [ $DO_NON_ITERATIVE -eq 0  -a  $DO_ITERATIVE -eq 0  -o  $MULTIDOMAINBENCHMARK_TRAINING -eq 0  -a  $MULTIDOMAINBENCHMARK_TEST -eq 0 ]; then
    echo "ERROR: Specify at least:" >&2
    echo "       -noniterative  or  -iterative" >&2
    echo "       and" >&2
    echo "       -training  or  -test" >&2
    echo "" >&2
    echo "" >&2
    echo "$usageStr" >&2
    exit 1
fi

calculateTaps=1
E_VALUE_THRESHOLD="1000"
#calculateTiming=0
calculateTiming=1

apps=""
appsIterations=""  # PSI-GLOBAL apps that have checkpoint PSSM for each initial search iteration
db=""


# check if MDB_RUNS_BASE_DIR is not an environmental variable
if [ -z "$MDB_RUNS_BASE_DIR" ]; then
    MDB_RUNS_BASE_DIR="MultiDomainBenchmarkRuns"
fi

# If MDB_DB_BASE_DIR is not an environmental variable, then set it to the currect directory
if [ -z "$MDB_DB_BASE_DIR" ]; then
    MDB_DB_BASE_DIR="."
fi

RELEVANCE_FILENAME="$MDB_DB_BASE_DIR/relevanceInfo.tab"



runJobsManagerAggregate(){
    if [ "$executeStr" == "-qsub" ]; then
	executeStr="-aggregate"
	runJobsManager "$@"
	executeStr="-qsub"
    else
	runJobsManager "$@"
    fi
} # END runJobsManagerAggregate()


runJobsManager(){

    if [ -z "$RELEVANCE_FILENAME" ]; then
	echo "ERROR: Relevance filename NOT SET!!!" >&2
	exit 1
    fi


    suffix=""
    flags=""
    if [ $# -lt 1 ] ; then
	echo "ERROR: runJobsManager() needs arg" >&2
	exit 1
    fi
    
    app=$1
    if [ $# -ge 2 ] ; then
	suffix=$2
	if [ $# -ge 3 ] ; then
	    flags=$3
	    if [ $# -ge 4 ] ; then
		flagsFinal=$4  # flags for the final iteration
	    fi
	fi
    fi
    
    echo "";
    echo "*****************************************************************************";
    echo "";
    #echo "suffix: $suffix";


    baseDirStr=""
    
    allDoneFile="$MDB_RUNS_BASE_DIR/$db/$app$suffix.allDone"
    if [ ! -f "$allDoneFile" ]; then
	flagsFileArg=""
	if [ -n "$flags" ]; then
	    flagsFile="$PWD/flagsFile_${app}_$suffix.txt"
	    flagsFileArg=" -flagsFile=$flagsFile "
	    if [ ! -s "$flagsFile" ]; then
		echo "$flags" > $flagsFile
	    else
		echo "$flags" > $flagsFile.tmp
		diffLine=`diff -q $flagsFile $flagsFile.tmp`
		if [ -n "$diffLine" ]; then
		    echo "" >&2
		    echo "diff -q $flagsFile $flagsFile.tmp" >&2
		    echo "" >&2
		    exit 1
		fi
		rm $flagsFile.tmp
	    fi
	fi
	flagsFinalFileArg=""
	if [ -n "$flagsFinal" ]; then
	    flagsFinalFile="$PWD/flagsFinalFile_${app}_$suffix.txt"
	    flagsFinalFileArg=" -flagsFinalFile=$flagsFinalFile "
	    if [ ! -s "$flagsFinalFile" ]; then
		echo "$flagsFinal" > $flagsFinalFile
	    else
		echo "$flagsFinal" > $flagsFinalFile.tmp
		diffLine=`diff -q $flagsFinalFile $flagsFinalFile.tmp`
		if [ -n "$diffLine" ]; then
		    echo "" >&2
		    echo "diff -q $flagsFinalFile $flagsFinalFile.tmp" >&2
		    echo "" >&2
		    exit 1
		fi
		rm $flagsFinalFile.tmp
	    fi
	fi

	cmd="./jobsManager.sh -$db -$app -suffix=$suffix $flagsFileArg $flagsFinalFileArg $countStr $executeStr $baseDirStr $pssmOnlyStr $iterationsDbStr $EVALUE_STR $verbosityStr "
	echo $cmd; $cmd
    fi

    
    if [ -f "$allDoneFile" ]; then
	VERBOSE "allDoneFile: \"$allDoneFile\" exists"
	if [ $calculateTaps -eq 1 ]; then
	    queryType="taxon"
	    cmd="./calculateTapKs.sh $MDB_RUNS_BASE_DIR/$db/$app$suffix $app$suffix $RELEVANCE_FILENAME $E_VALUE_THRESHOLD $queryType $db"
	    VERBOSE "$cmd"
	    $cmd

	    if [ "$app" == "psiblast"  -o  "$app" == "psisemiglobal" ]; then
		startIteration=1
		if [ "$app" == "psiblast" ]; then
		    startIteration=2  # psiblast -query ... does not have a PSSM for iteration 1
		fi
	    	cmd="./calculateTap-iterations.sh $MDB_RUNS_BASE_DIR/$db/$app$suffix $app$suffix $RELEVANCE_FILENAME $E_VALUE_THRESHOLD $queryType $db  $startIteration"
	    	VERBOSE "$cmd"
	    	$cmd
	    	appsIterations="$appsIterations $app$suffix"
	    fi
	    
	    apps="$apps $app$suffix"
	fi
	if [ $calculateTiming -eq 1 ]; then
	    cmd="./calculateTiming.sh $MDB_RUNS_BASE_DIR/$db/$app$suffix $app$suffix"
	    VERBOSE "$cmd"
	    $cmd
	fi
    fi
} # END runJobsManager


stats(){
    echo "Collecting TAP-ks and timing values . . ."

    descriptor="$1"
    
    timingOutFilename="aggregateTiming$descriptor.csv"
    if [ -s "$timingOutFilename" ]; then
	mv -f "$timingOutFilename" "$timingOutFilename.bak"
    fi

    # remove duplicate apps
    apps=`for app in $apps; do echo "$app"; done | sort -n | uniq`

    str="DB"
    for app in $apps; do
	str="$str,$app"
    done
    echo $str > "$timingOutFilename"

    tapKsOutFilename="aggregateTapKs$descriptor.csv"
    if [ -s "$tapKsOutFilename" ]; then
	mv -f "$tapKsOutFilename" "$tapKsOutFilename.bak"
    fi
    
    str="DB,TAP"
    for app in $apps; do
	str="$str,$app"
    done
    echo $str > "$tapKsOutFilename"

    ks="1 3 5 20"
    for db in $dbs; do
	for k in $ks; do
	    str="$db,TAP-$k"
	    for app in $apps; do
		tapKFilename="$MDB_RUNS_BASE_DIR/$db/$app.tap$k"
		if [ ! -s "$tapKFilename" ]; then
		    str="$str,0"
		else
		    str="$str,"`cat $tapKFilename`
		    echo "$tapKFilename: "`cat $tapKFilename`
		fi
	    done
	    echo $str >> "$tapKsOutFilename"
	done

        # timing
	str="$db"
	for app in $apps; do
	    timingFilename="$MDB_RUNS_BASE_DIR/$db/$app.timing.out"
	    if [ ! -s "$timingFilename" ]; then
		str="$str,0"
	    else
		str="$str,"`cat $timingFilename`
	    fi
	done
	echo $str >> "$timingOutFilename"
    done
} # END stats()



stats-iterations(){
    echo "Collecting TAP-ks and timing values (iterative) . . ."
    echo "appsIterations: $appsIterations"
        
    descriptor="$1"
    
    tapKsOutFilename="aggregateTapKs-iterations$descriptor.csv"
    if [ -s "$tapKsOutFilename" ]; then
	mv -f "$tapKsOutFilename" "$tapKsOutFilename.bak"
    fi
    
    str="DB,TAP-k"
    k=1
    for db in $dbs; do
	for app in $appsIterations; do
	    for iteration in 1 2 3 4 5; do
		tapKFilename="$MDB_RUNS_BASE_DIR/$db/$app.tap${k}_$iteration"
		if [ -s "$tapKFilename" ]; then
		    str="$str,$app (iter $iteration)"
		else
		    str="$str,"
		fi
	    done
	done
	break;  # just do 1 DB
    done
    
    echo $str > "$tapKsOutFilename"

    ks="1 3 5 20"
    for db in $dbs; do
	for k in $ks; do
	    str="$db,TAP-$k"
	    for app in $appsIterations; do
		for iteration in 1 2 3 4 5; do
		    tapKFilename="$MDB_RUNS_BASE_DIR/$db/$app.tap${k}_$iteration"
		    if [ -s "$tapKFilename" ]; then
			str="$str,"`cat $tapKFilename`
			echo "$tapKFilename: "`cat $tapKFilename`
		    else
			str="$str,"
		    fi
		done
	    done
	    echo $str >> "$tapKsOutFilename"
	done
    done




    tapKsOutFilename="aggregateTapKs-iterations2$descriptor.csv"
    if [ -s "$tapKsOutFilename" ]; then
	mv -f "$tapKsOutFilename" "$tapKsOutFilename.bak"
    fi
    
    ks="1 3 5 20"
    for db in $dbs; do
        # header
	str="$db"
	for k in $ks; do
	    str="$str,TAP-$k"
	    for iteration in 1 2 3 4 5; do
		str="$str,"
	    done
	done
	echo $str >> "$tapKsOutFilename"

	str="Iteration"
	for k in $ks; do
	    for iteration in 1 2 3 4 5; do
		str="$str,$iteration"
	    done
	    str="$str,"
	done
	echo $str >> "$tapKsOutFilename"

	# data
	for app in $appsIterations; do
	    str="$app"
	    for k in $ks; do
		for iteration in 1 2 3 4 5; do
		    tapKFilename="$MDB_RUNS_BASE_DIR/$db/$app.tap${k}_$iteration"
		    if [ -s "$tapKFilename" ]; then
			str="$str,"`cat $tapKFilename`
		    else
			str="$str,"
		    fi
		done
		str="$str,"
	    done
	    echo "$str" >> "$tapKsOutFilename"
	done
    done
} # END stats-iterations()



function VERBOSE(){
    if [ $VERBOSE -eq 1 ]; then
	echo "$1"
    fi
} # END VERBOSE


function DEBUG(){
    if [ $DEBUG -eq 1 ]; then
	echo "$1"
    fi
} # END DEBUG

# END subroutines ##################################################################



if [ $MULTIDOMAINBENCHMARK_TRAINING -eq 1 ]; then
    db="MultiDomainBenchmark-training"
    dbs="$dbs $db"

    DEBUG "$db";

    if [ $calculateTaps -eq 1 ]; then

	if [ ! -s "$RELEVANCE_FILENAME" ]; then
	    echo "ALERT: Relevance file: \"$RELEVANCE_FILENAME\" does not exist!!! It's needed to populate empty retrieval files." >&2
	    exit 1
	fi
    fi
    
    if [ $DO_NON_ITERATIVE -eq 1 ]; then
        #
        # non_iterative
        #
	
	#runJobsManagerAggregate "psiblast-non_iterative" "" "$BLAST_COMMON_ARGS"
	runJobsManager "psiblast-non_iterative" "" "$BLAST_COMMON_ARGS"
	#runJobsManager "psiblast-non_iterative" "-defaults"

	#len="35"
	#for len in 20 25 30 35 40 45; do
	for len in 20 30 40; do
	    flag="-target_block_length $len"
	    flags=" $flag  $BLAST_COMMON_ARGS"
	    #runJobsManagerAggregate
	    runJobsManager "psisemiglobal-non_iterative" "-blast_pssm-targetBlockLen$len" "$flags"
	done

	# len=30
	# flag="-target_block_length $len"
	# flags=" $flag " ##### $BLAST_COMMON_ARGS"
	# runJobsManager "psisemiglobal-non_iterative" "-defaults-blast_pssm-targetBlockLen$len" "$flags"

	#for len in 20 25 30 35 40 45; do
	for len in 20 30 40; do
	    flag="-target_block_length $len"
	    flags=" -global_pssm  $flag  $BLAST_COMMON_ARGS"
	    #runJobsManagerAggregate
	    runJobsManager "psisemiglobal-non_iterative" "-global_pssm-targetBlockLen$len" "$flags"
	done

	# len=30
	# flag="-target_block_length $len"
	# flags=" -global_pssm  $flag " ##### $BLAST_COMMON_ARGS"
	# runJobsManager "psisemiglobal-non_iterative" "-defaults-global_pssm-targetBlockLen$len" "$flags"

	#runJobsManager "global"
	
    fi  # END if [ $DO_NON_ITERATIVE -eq 1 ]


    if [ $DO_ITERATIVE -eq 1 ]; then
        #
        # iterative
        #
 	runJobsManager "psiblast"  ""  ""  "$BLAST_COMMON_ARGS"
	#runJobsManager "psiblast"  "-defaults"
	
	#len="35"
	#for len in 20 25 30 35 40 45; do
	for len in 20 30 40; do
	    setFlags=" -target_block_length $len"
	    runJobsManager "psisemiglobal" "-iterative-blast_pssm-targetBlockLen$len" "$setFlags" "$setFlags  $BLAST_COMMON_ARGS"

	    setFlags=" -target_block_length $len  -update_block_indices"
	    runJobsManager "psisemiglobal" "-iterative-blast_pssm-targetBlockLen$len-updateBlockIndices" "$setFlags" "$setFlags  $BLAST_COMMON_ARGS"
	done

	# len=30
	# setFlags=" -target_block_length $len "
	# runJobsManager "psisemiglobal" "-iterative-defaults-blast_pssm-targetBlockLen$len" "$setFlags" "$setFlags" #  $BLAST_COMMON_ARGS"

	# setFlags=" -target_block_length $len  -update_block_indices"
	# runJobsManager "psisemiglobal" "-iterative-defaults-blast_pssm-targetBlockLen$len-updateBlockIndices" "$setFlags" "$setFlags" #  $BLAST_COMMON_ARGS"


	#for len in 20 25 30 35 40 45; do
	for len in 20 30 40; do
	    setFlags=" -global_pssm  -target_block_length $len"
	    runJobsManager "psisemiglobal" "-iterative-global_pssm-targetBlockLen$len" "$setFlags" "$setFlags  $BLAST_COMMON_ARGS"

	    setFlags=" -global_pssm  -target_block_length $len  -update_block_indices"
	    runJobsManager "psisemiglobal" "-iterative-global_pssm-targetBlockLen$len-updateBlockIndices" "$setFlags" "$setFlags  $BLAST_COMMON_ARGS"
	done

	# len=30
	# setFlags=" -global_pssm  -target_block_length $len "
	# runJobsManager "psisemiglobal" "-iterative-defaults-global_pssm-targetBlockLen$len" "$setFlags" "$setFlags" #  $BLAST_COMMON_ARGS"

	# setFlags=" -global_pssm  -target_block_length $len  -update_block_indices"
	# runJobsManager "psisemiglobal" "-iterative-defaults-global_pssm-targetBlockLen$len-updateBlockIndices" "$setFlags" "$setFlags" #  $BLAST_COMMON_ARGS"


        # runJobsManager "jackhmmer"
    fi  # END if [ $DO_ITERATIVE -eq 1 ]
fi # END [ $MULTIDOMAINBENCHMARK_TRAINING -eq 1 ]


if [ $MULTIDOMAINBENCHMARK_TEST -eq 1 ]; then
    db="MultiDomainBenchmark-test"
    dbs="$dbs $db"

    DEBUG "$db";

    if [ $calculateTaps -eq 1 ]; then
	if [ ! -s "$RELEVANCE_FILENAME" ]; then
	    echo "ALERT: Relevance filename: \"$RELEVANCE_FILENAME\" does not exist!!! It's needed to populate empty retrieval files." >&2
	    exit 1
	fi
    fi
 
    if [ $DO_NON_ITERATIVE -eq 1 ]; then
        #
        # non_iterative
        #
	
	#runJobsManagerAggregate "psiblast-non_iterative" "" "$BLAST_COMMON_ARGS"
	runJobsManager "psiblast-non_iterative" "" "$BLAST_COMMON_ARGS"
	#runJobsManager "psiblast-non_iterative" "-defaults"  # uncomment this line and comment out the line above if not worried about ROC values 
	
        len="30"
	flag="-target_block_length $len"
	# flags=" -global_pssm  $flag " # $BLAST_COMMON_ARGS"
	# runJobsManager "psisemiglobal-non_iterative" "-defaults-global_pssm-targetBlockLen$len" "$flags"
	
	flags=" -global_pssm  $flag  $BLAST_COMMON_ARGS"
	runJobsManager "psisemiglobal-non_iterative" "-global_pssm-targetBlockLen$len" "$flags"

	# HMMER
	runJobsManager "phmmer"
    fi  # END if [ $DO_NON_ITERATIVE -eq 1 ]


    if [ $DO_ITERATIVE -eq 1 ]; then
        #
        # iterative
        #
	#runJobsManager "global"

 	runJobsManager "psiblast"  ""  ""  "$BLAST_COMMON_ARGS"
	#runJobsManager "psiblast"  "-defaults"


        len="30"
	flag="-target_block_length $len"
	# flags=" -global_pssm  $flag " # $BLAST_COMMON_ARGS"
	# runJobsManager "psisemiglobal" "-iterative-defaults-global_pssm-targetBlockLen$len" "$flags"

	flags=" -global_pssm  $flag"
	runJobsManager "psisemiglobal"  "-iterative-global_pssm-targetBlockLen$len"  "$flags"  "$flags  $BLAST_COMMON_ARGS"
	
        runJobsManager "jackhmmer"
    fi  # END if [ $DO_ITERATIVE -eq 1 ]
fi # if [ $MULTIDOMAINBENCHMARK_TEST -eq 1 ]



# descriptor for filenames
dbStr=""
for db in $dbs; do
    if [ -n "$dbStr" ]; then
	    # put "_"s inbetween entries
	dbStr="${dbStr}__"
    fi
    dbStr="$dbStr$db"	    
done
if [ -n "$dbStr" ]; then
    dbStr="-$dbStr"
fi

descriptor=""
if [ $DO_NON_ITERATIVE -eq 1 ]; then
    descriptor="non_iterative"
fi
if [ $DO_ITERATIVE -eq 1 ]; then
    if [ -n "$descriptor" ]; then
	    # put "_"s inbetween entries
	descriptor="${descriptor}__"
    fi
    descriptor="${descriptor}iterative"
fi
if [ -n "$descriptor" ]; then
    descriptor="-$descriptor"
fi
descriptor="$dbStr$descriptor"
    

stats $descriptor
if [ -n "$appsIterations" ]; then
    stats-iterations $descriptor
fi


exit

