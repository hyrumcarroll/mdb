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


# Note: UPPERCASE variables are global consts or commandline parameters

# make a separate .spouge file for ROC_n analysis?
#ROCN_SPOUGE_FILES=1
ROCN_SPOUGE_FILES=0


WALLTIME="36000"  # number of seconds for submitted jobs (36000s = 10 hours)
AGGREGATE_JOBS_NUMBER=50  # number of jobs to lump together in a single qsub 
#AGGREGATE_JOBS_NUMBER=150  # number of jobs to lump together in a single qsub 
AGGREGATE_JOBS_COUNT=$AGGREGATE_JOBS_NUMBER
AGGREGATE_JOBS_SUBMITTED_COUNTER=0

scriptDir=${0%/*}
if [ "$scriptDir" == "$0" ]; then
    scriptDir="."
fi
scriptLocalFile=${0##*/}
LOCK_FILE="$scriptDir/.lock.$scriptLocalFile"

# check if MDB_RUNS_BASE_DIR is not an environmental variable
if [ -z "$MDB_RUNS_BASE_DIR" ]; then
    MDB_RUNS_BASE_DIR="MultiDomainBenchmarkRuns"
fi
if [ ! -d $MDB_RUNS_BASE_DIR ]; then
    cmd="mkdir -p $MDB_RUNS_BASE_DIR"
    echo $cmd
    $cmd
fi

# If MDB_DB_BASE_DIR is not an environmental variable, then set it to the currect directory
if [ -z "$MDB_DB_BASE_DIR" ]; then
    MDB_DB_BASE_DIR="."
else
    if [ ! -d "$MDB_DB_BASE_DIR" ]; then
	echo "ERROR: MultiDomainBenchmark database base directory \"$MDB_DB_BASE_DIR\" does not exist!" >&2
	exit 1
    fi
fi

# Set MDB_BENCHMARKING_SCRIPTS_DIR to the currently directory unless it is an environmental variable
if [ -z "$MDB_BENCHMARKING_SCRIPTS_DIR" ]; then
    MDB_BENCHMARKING_SCRIPTS_DIR="."
else
    if [ ! -d "$MDB_BENCHMARKING_SCRIPTS_DIR" ]; then
	echo "ERROR: MultiDomainBenchmark benchmarking scripts directory \"$MDB_BENCHMARKING_SCRIPTS_DIR\" does not exist!" >&2
	exit 1
    fi
fi


SLEEP_COUNTER=0

VERBOSE=0
DEBUG=0
argsStr=""
QSUB=0
SH=0
AGGREGATE_JOBS=0
MAKE_SCRIPTS=0
SCORE_ONLY=0
RESUBMIT=0
TO_QUEUE_COUNT=1    # submit one job
TOTAL_QUEUE_COUNT=0 # submit jobs until there is 1 job in the queue
SUFFIX=""
FLAGS_FILE=""
FLAGSFINAL_FILE=""
E_VALUE_THRESHOLD=""
PSSM_ONLY=0
ITERATIVE=1
ITERATIONS_DB=""
ITERATIONS_DB_STR=""

FINISHED_LINE="Finished - ALL DONE"


usage(){
    echo "Usage: $0 [parameters]
    	# applications:
	-psisemiglobal
	-psiblast 
#	-blastp
	-psiblast-non_iterative
	-global   
	-globalNrTiming
	-phmmer
	-jackhmmer
	-allapps  
    
        # benchmark DBs:
	-MultiDomainBenchmark-training
	-MultiDomainBenchmark-test
	-alldbs

	-all   
    
        -scripts | -script
	-qsub  
	-sh   
	-score 

	-queries=* | -query=* | -q=*
	
	-queuecount=* | -queue=*
        -totalcount=* | -total=*

	-suffix=*

	-flagsFile=<file> | -flags=<file> | -flag=<file> | -f=<file>

	-flagsFinalFile=<file> | -flagsFinal=<file>

	-evalue=<number (evalue threshold)>

        -iterationsDb <file (full path filename to database)>

	-basedir <file (MultiDomainBenchmark runs base directory)>

	-pssmOnly

	-verbose | -v
	-debug | -d

	-test | -t
	" >&2
}

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

submitAggregateJobs(){
    tag="$SUFFIX"
    if [ -z "$tag" ]; then
	tag="$$"
    fi
    aggregateSubmitScriptFile="aggregate${tag}-$AGGREGATE_JOBS_SUBMITTED_COUNTER.sh"
    if [ $DEBUG -eq 0 ]; then   rm -f $aggregateSubmitScriptFile;   fi
    AGGREGATE_JOBS_SUBMITTED_COUNTER=$(( $AGGREGATE_JOBS_SUBMITTED_COUNTER + 1 ))
    echo "#!/bin/bash
	
#$ -N aggregate
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -j n
#$ -P unified
#$ -l h_rt=$WALLTIME
#$ -l h_vmem=8.1G
#$ -l mem_free=6G
#$ -m n

source ~/.bashrc
ulimit -v $(( 8 * 1024 * 1024 ))

$AGGREGATE_JOBS_STR
" > $aggregateSubmitScriptFile
    qsub $aggregateSubmitScriptFile
    #rm $aggregateSubmitScriptFile
    AGGREGATE_JOBS_STR=""
    AGGREGATE_JOBS_COUNT=$AGGREGATE_JOBS_NUMBER
    sleep 1
} # END submitAggregateJobs()

submitJob(){
    allDone=0
    if [ $remainingJobsToQueue -eq 0  ]; then
	#return
	cleanUpAndExit 0
    fi
    
    if [ $MAKE_SCRIPTS -eq 1 ]; then
	submitFile=$submitDir"/"$jobName".sh"
	
	rm -f $submitFile

	submitScript="#!/bin/bash
	
#$ -N $jobName
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -j n
## -o $outFile.job
## -e $errFile.job
#$ -P unified
#$ -l h_rt=$WALLTIME
#$ -l h_vmem=8.1G
#$ -l mem_free=6G
# use reserve_mem only for large jobs b/c it will block other users from getting access to that node
#$ -l reserve_mem=8G
## -m a
#$ -m n


source ~/.bashrc
ulimit -v $(( 8 * 1024 * 1024 ))

jobName=\"$jobName\"
queryFile=\"$queryFile\"
spougeFile=\"$spougeFile\"
outFile=\"$outFile\"
errFile=\"$errFile\"
iterative=$ITERATIVE
hitsFile=\"$finalHitsFile\"
outputFile=\"$outputFile\"
if [ \$iterative -eq 1 ]; then 
    pssmFile=\"$pssmFile\"
    if [ -n \"\$pssmFile\" ]; then
        pssmStr=\" -pssm=\$pssmFile \"
        pssmOnly=$PSSM_ONLY
        if [ \$pssmOnly -eq 1 ]; then
            pssmStr=\" \$pssmStr -pssmOnly \"
        fi
    fi
    iterativeStr=\" \$pssmStr  $ITERATIONS_DB_STR \"
fi

cd $submitDir

if [ ! -d \"\$jobName\" ]; then
    cmd=\"mkdir \$jobName\"
    #echo \"\$cmd\"
    \$cmd
fi
cd \$jobName

for filename in \$outFile \$errFile; do
    if [ -f \$filename.gz ]; then
        gunzip \$filename.gz
    fi
    echo \"-----------------------------------------------------------------------------\"\`date\` >> \$filename
done    

cmd=\"$appWrapper -$app -query=\$queryFile  -spouge=\$spougeFile -db=$finalDatabase  -rel=$relevanceFile  $argsStr  $flagsFileStr  $flagsFinalFileStr   -hits=\$hitsFile -out=\$outputFile  $evalueStr  \$iterativeStr \" 

echo \"\$cmd\" >> \$outFile 2>> \$errFile
if [ $TEST -eq 1 ]; then exit; fi
\$cmd >> \$outFile 2>> \$errFile


rc=\$?
if [ \$rc != 0 ]; then
    echo \":\$cmd: FAILED (rc = \$rc)\" >> \$errFile
    exit 1
fi

if [ -n \"\$pssmFile\" ]; then
    if [ \$pssmOnly -eq 1 ]; then
        if [ ! -s \$pssmFile ]; then
            echo \"PSSM file, \$pssmFile, does not exist!\" >> \$errFile
            exit 1
        fi
    elif [ ! -s \$pssmFile.gz ]; then
            echo \"PSSM file, \$pssmFile.gz, does not exist!\" >> \$errFile
            exit 1
    fi
fi

if [ ! -s \$spougeFile ]; then
    echo \"Spouge file, \$spougeFile, does not exist!\" >> \$errFile
    exit 1
fi

if [ ! -s \$errFile ]; then
    rm \$errFile
fi

cd $submitDir

(echo \"\"; echo \"$FINISHED_LINE\") >> \$outFile
"

	if [ $remainingJobsToQueue -eq 1  -a  $QSUB -eq 1  -a  $RESUBMIT -ge 1 ]; then
	    submitScript="$submitScript
$MDB_BENCHMARKING_SCRIPTS_DIR/jobsManager.sh -qsub -total=$RESUBMIT  $resubmitStr"
	    RESUBMIT=0
	fi

	submitScript="$submitScript
exit 0"
	echo "$submitScript" > $submitFile

	if [ $QSUB -eq 1 ]; then
            #cmd="qsub $submitFile -o /dev/null -e /dev/null"
	    cmd="qsub $submitFile -o $outFile.job -e $errFile.job"
	    if [ $VERBOSE -eq 1 ]; then echo "$cmd"; fi
	    $cmd

	    SLEEP_COUNTER=$(( $SLEEP_COUNTER + 1 ))
	    if [ $SLEEP_COUNTER -eq 5 ]; then
		SLEEP_COUNTER=0
		cmd="sleep 1"
	        #if [ $VERBOSE -eq 1 ]; then
	    	    echo "$cmd";
	        #fi
		$cmd
	    fi
	elif [ $SH -eq 1 ]; then
	    cmd="sh $submitFile"
	    if [ $VERBOSE -eq 1 ]; then echo "$cmd"; fi
	    $cmd
	elif [ $AGGREGATE_JOBS -eq 1 ]; then
	    AGGREGATE_JOBS_STR="$AGGREGATE_JOBS_STR 
sh $submitFile"
	    AGGREGATE_JOBS_COUNT=$(( $AGGREGATE_JOBS_COUNT - 1 ))
	    if [ $AGGREGATE_JOBS_COUNT -eq 0 ]; then
		submitAggregateJobs
	    fi
	fi
    fi # if( MAKE_SCRIPTS)
    remainingJobsToQueue=$(( $remainingJobsToQueue - 1))	
    #if [ $VERBOSE -eq 1 ]; then echo "remainingJobsToQueue: $remainingJobsToQueue"; fi
    if [ $remainingJobsToQueue -eq 0  ]; then
	#return
	cleanUpAndExit 0
    fi
} # END submitJob()


# @return 0 if successful, 1 if not
function setLock(){
    numTries=100
    tryCounter=0
    result=-1
    
    while [ $tryCounter -lt $numTries ]; do
	if [ $DEBUG -eq 1 ]; then echo "Try: $tryCounter ($$)"; fi
	if set -o noclobber; echo "$$" > "$LOCK_FILE" 2> /dev/null ; then
            # writing to lock file is successful
            result=0
	    echo "Wrote to LOCK_FILE \"$LOCK_FILE\" ($$)"
	    break
	fi

	tryCounter=$(( $tryCounter + 1 ))
	sleep 3
    done

    if [ $result -ne 0 ]; then
        # writing to lock file yields errors;
        # Does lock file hold a valid pid?
        if test -n "$(cat "$LOCK_FILE")" && ps p $(cat "$LOCK_FILE") > /dev/null; then
            # Yes, lock file holds a valid pid;
            if [ $DEBUG -eq 1 ]; then echo "LOCK_FILE already exists ($$)"; fi
            result=1
            echo "ERROR: LOCK_FILE, $LOCK_FILE, set by existing process (pid = $(cat $LOCK_FILE));"  >&2
	    exit 1
        else
            # No, lock file holds a non-valid pid
            # Indicate error
            echo "ERROR: LOCK_FILE, $LOCK_FILE, set by non-existing process (pid = $(cat $LOCK_FILE)); removing lock ($$)"  >&2
	    unsetLock $LOCK_FILE
            cleanUpAndExit 1
        fi
    fi
    return $result
} # setlock()

function unsetLock(){
    # Remove lockfile
    #if [ $DEBUG -eq 1 ]; then echo "Removing LOCK_FILE \"$LOCK_FILE\" ($$)"; fi
    echo "Removing LOCK_FILE \"$LOCK_FILE\" ($$)"
    rm -f $LOCK_FILE
} # unsetLock()

cleanUpAndExit(){
    exitCode=$1

    if [ $AGGREGATE_JOBS -eq 1  -a  $AGGREGATE_JOBS_COUNT -gt 0  -a $AGGREGATE_JOBS_COUNT -ne $AGGREGATE_JOBS_NUMBER ]; then
	submitAggregateJobs
    fi

	    
    /bin/rm $jobsFile
    unsetLock
    
    if [ -z "$exitCode" ]; then
	exitCode=0
    fi

    exit $exitCode
} # cleanUp()



# Check for file indicating no more scheduling
NO_MORE_SCHEDULING_FILE="$HOME/noMoreScheduling.txt"
if [ -f "$NO_MORE_SCHEDULING_FILE" ]; then
    echo "NO_MORE_SCHEDULING" >&2
    exit
fi

resubmitStr=""

TEST=0
DBS=""

ALL_APPS="psisemiglobal psiblast global globalNrTiming psisemiglobal-non_iterative psiblast-non_iterative phmmer jackhmmer"
ALL_DBS="MultiDomainBenchmark"

for arg in "$@" ; do
    #echo "arg = $arg"
    case "$arg" in
	# applications:
	-psisemiglobal   )  APPS="$APPS psisemiglobal" ;;
	-psiblast    )  APPS="$APPS psiblast" ;;
	-psiblast-non_iterative )  APPS="$APPS psiblast-non_iterative"; ITERATIVE=0 ;;
	-global      )  APPS="$APPS global"; ITERATIVE=0 ;;
	-globalNrTiming      )  APPS="$APPS globalNrTiming"; ITERATIVE=0 ;;
	-psisemiglobal-non_iterative   )  APPS="$APPS psisemiglobal-non_iterative"; ITERATIVE=0 ;;
	-phmmer      )  APPS="$APPS phmmer" ; ITERATIVE=0 ;;
	-jackhmmer   )  APPS="$APPS jackhmmer" ;;
	#-blastp       )  APPS="$APPS blastp"; ITERATIVE=0 ;;
	-allapps     )  APPS="$ALL_APPS" ;;
    
        # benchmark DBs:
	-MultiDomainBenchmark-training | -MDB-training | -mdb-training )  DBS="$DBS MultiDomainBenchmark-training" ;;
	-MultiDomainBenchmark-test | -MDB-test | -mdb-test )  DBS="$DBS MultiDomainBenchmark-test" ;;

#	-nr          )  DBS="$DBS nr" ;;
	
	-iterationsDb=* )  ITERATIONS_DB_STR="$arg" ;;

	-alldbs      )  DBS="$ALL_DBS" ;;

	-all         )  APPS="$ALL_APPS"; DBS="$ALL_DBS" ;;

    
        -scripts     )  MAKE_SCRIPTS=1 ;;
        -script      )  MAKE_SCRIPTS=1 ;;
	-qsub        )  QSUB=1; MAKE_SCRIPTS=1 ;;
	-sh          )  SH=1;   MAKE_SCRIPTS=1 ;;
	-aggregate   )  AGGREGATE_JOBS=1;   MAKE_SCRIPTS=1 ;;
	-score       )  SCORE_ONLY=1; resubmitStr="$resubmitStr $arg" ;;
	-resubmit    )  RESUBMIT=1;   resubmitStr="$resubmitStr $arg" ;;

	-queries=*   )  QUERIES="`echo $arg | sed -e 's/-queries=//'`" ;;
	-query=*     )  QUERIES="`echo $arg | sed -e 's/-query=//'`" ;;
	-q=*         )  QUERIES="`echo $arg | sed -e 's/-q=//'`" ;;
	
# 	-count=*     )  COUNT="`echo $arg | sed -e 's/-count=//'`" ;;
# 	-c=*         )  COUNT="`echo $arg | sed -e 's/-c=//'`" ;;
	-queue=*     )  TO_QUEUE_COUNT="`echo $arg | sed -e 's/-queue=//'`" ;;
	-queuecount=*)  TO_QUEUE_COUNT="`echo $arg | sed -e 's/-queuecount=//'`" ;;
	-total=*     )  TOTAL_QUEUE_COUNT="`echo $arg | sed -e 's/-total=//'`" ;;
	-totalcount=*)  TOTAL_QUEUE_COUNT="`echo $arg | sed -e 's/-totalcount=//'`" ;;

	-suffix=*    )  SUFFIX="`echo $arg | sed -e 's/-suffix=//'`"; resubmitStr="$resubmitStr $arg" ;;
	
	-flagsFile=* )  FLAGS_FILE="`echo $arg | sed -e 's/-flagsFile=//'`";  resubmitStr="$resubmitStr $arg" ;;
	-flagFile=*  )  FLAGS_FILE="`echo $arg | sed -e 's/-flagFile=//'`";   resubmitStr="$resubmitStr $arg" ;;
	-flags=*     )  FLAGS_FILE="`echo $arg | sed -e 's/-flags=//'`";      resubmitStr="$resubmitStr $arg" ;;
	-flag=*      )  FLAGS_FILE="`echo $arg | sed -e 's/-flag=//'`";       resubmitStr="$resubmitStr $arg" ;;
	-f=*         )  FLAGS_FILE="`echo $arg | sed -e 's/-f=//'`";          resubmitStr="$resubmitStr $arg" ;;
	
	-flagsFinalFile=* )  FLAGSFINAL_FILE="`echo $arg | sed -e 's/-flagsFinalFile=//'`";  resubmitStr="$resubmitStr $arg" ;;
	-flagsFinal=*     )  FLAGSFINAL_FILE="`echo $arg | sed -e 's/-flagsFinal=//'`";      resubmitStr="$resubmitStr $arg" ;;
	
	-evalue=*    )  E_VALUE_THRESHOLD="`echo $arg | sed -e 's/-evalue=//'`"; resubmitStr="$resubmitStr $arg" ;;
	-eValue=*    )  E_VALUE_THRESHOLD="`echo $arg | sed -e 's/-eValue=//'`"; resubmitStr="$resubmitStr $arg" ;;
	-e=*         )  E_VALUE_THRESHOLD="`echo $arg | sed -e 's/-e=//'`"; resubmitStr="$resubmitStr $arg" ;;

	-basedir=*   )  MDB_RUNS_BASE_DIR="`echo $arg | sed -e 's/-basedir=//'`" ;;
	
        -pssmOnly    )  PSSM_ONLY=1 ;;
			
	-verbose     )  VERBOSE=1 ;;
	-v           )  VERBOSE=1 ;;
	-debug       )  DEBUG=1 ;; 
	-d           )  DEBUG=1 ;;

	-test        )  TEST=1 ;;
	-t           )  TEST=1 ;;

	-help        )  usage; exit 0 ;;
	--help       )  usage; exit 0 ;;
			
	*            )  echo "" >&2; echo "ALERT: Skipping unrecgonized parameter: \"$arg\"" >&2; echo "" >&2; ;;
  esac
done

#DEBUG_FILE_STR=""   # appended onto files/directories

if [ $DEBUG -eq 1 ]; then
    VERBOSE=1
    argsStr="$argsStr -d"
    #DEBUG_FILE_STR=".debug"
fi

if [ $VERBOSE -eq 1 ]; then
    argsStr="$argsStr -v"
fi

if [ $RESUBMIT -eq 1 ]; then
    RESUBMIT=$TOTAL_QUEUE_COUNT
fi

resubmitStr="$resubmitStr $argsStr"

for app in $APPS; do
    resubmitStr="$resubmitStr -$app"
done

for db in $DBS; do
    resubmitStr="$resubmitStr -$db"
done

if [ -z "$DBS" ]; then
#     #DBS="astral"
#     #DBS="cdd"
#     echo "NOTE: No database parameters set, using $DBS"
    usage
    exit 1
fi

flagsFileStr=""
if [ -n "$FLAGS_FILE" ]; then
    flagsFileStr="-flagsFile=$FLAGS_FILE"
fi

flagsFinalFileStr=""
if [ -n "$FLAGSFINAL_FILE" ]; then
    flagsFinalFileStr="-flagsFinalFile=$FLAGSFINAL_FILE"
fi

if [ ! -d $MDB_RUNS_BASE_DIR ]; then
    echo "ERROR: MDB_RUNS_BASE_DIR \"$MDB_RUNS_BASE_DIR\" does not exist" 1>&2
    exit 1
fi

if [ $VERBOSE -eq 1 ]; then
    echo ""
    echo "$0 args:"
    echo "  DBS:                           $DBS"
    echo "  APPS:                          $APPS"
    echo "  QUERIES(should be full path):   $QUERIES"
    echo "  TEST:                           $TEST (if set: print commandline, but don't execute runs)"
    echo "  MAKE_SCRIPTS:                   $MAKE_SCRIPTS"
    echo "  QSUB:                           $QSUB"
    echo "  SH:                             $SH"
    echo "  AGGREGATE_JOBS:                 $AGGREGATE_JOBS"
    #echo "  SCORE_ONLY:                     $SCORE_ONLY"
    echo "  TO_QUEUE_COUNT:                 $TO_QUEUE_COUNT"
    echo "  TOTAL_QUEUE_COUNT:              $TOTAL_QUEUE_COUNT"
    echo "  SUFFIX:                         $SUFFIX"
    echo "  FLAGS_FILE:                     $FLAGS_FILE"
    echo "  FLAGSFINAL_FILE:                $FLAGSFINAL_FILE"
    echo "  E_VALUE_THRESHOLD:              $E_VALUE_THRESHOLD"
    echo "  MDB_RUNS_BASE_DIR:              $MDB_RUNS_BASE_DIR"
    echo "  ITERATIVE:                      $ITERATIVE"
    echo "  ITERATIONS_DB_STR:              $ITERATIONS_DB_STR"
    echo "  Verbose:                        $VERBOSE"
    echo "  Debug:                          $DEBUG"
fi

evalueStr="" # for submit script
if [ $ITERATIVE -eq 0  -o  $PSSM_ONLY -ne 1 ]; then
    if [ -n "$E_VALUE_THRESHOLD" ]; then
	evalueStr=" -evalue=$E_VALUE_THRESHOLD "
    fi
fi
	        
# check for invalid query files
for queryFile in $QUERIES; do
    exitIfFileDoesNotExist $queryFile "ERROR: Invalid query file \"$queryFile\""
done

# if -total or -totalcount is set, then override -queue | -queuecount
if [ $TOTAL_QUEUE_COUNT -gt 0 ]; then
    TO_QUEUE_COUNT=0;
fi


#
# Get SGE jobs
#
setLock

jobsFile="$MDB_RUNS_BASE_DIR/.queuedAndRunningJobs.$$.xml"
qstat -u $USER -xml | grep JB_name > $jobsFile
perl -pi -e 's/^\s*\<JB_name>([^\<]+)\<\/JB_name\>/$1/' $jobsFile
if [ $DEBUG -eq 1 ]; then cmd="cat $jobsFile"; echo $cmd; $cmd;  fi

#
# calculate the number of jobs to queue based on TO_QUEUE_COUNT or TOTAL_QUEUE_COUNT and what's already queued
#
jobsInQueue=`cat $jobsFile | wc -l`
if [ $DEBUG -eq 1 ]; then echo "jobsInQueue: $jobsInQueue"; fi
if [ $TO_QUEUE_COUNT -gt 0 ]; then
    remainingJobsToQueue=$TO_QUEUE_COUNT
else
    if [ $QSUB -eq 1 ]; then # -o  $SH -eq 1
	totalJobsToQueue=$(( TOTAL_QUEUE_COUNT - jobsInQueue ))
	if [ $totalJobsToQueue -le 0 ]; then
	    echo "Requested up to $TOTAL_QUEUE_COUNT jobs to queue, but there's too many ($jobsInQueue) already in the queue"
	    cleanUpAndExit 0
	else
	    if [ $VERBOSE -eq 1 ]; then echo "totalJobsToQueue: $totalJobsToQueue"; fi
	fi
	remainingJobsToQueue=$totalJobsToQueue
    else
	# we're not submitting a job, so use the full TOTAL_QUEUE_COUNT
	remainingJobsToQueue=$TOTAL_QUEUE_COUNT
    fi
fi

    
for app in $APPS; do
    echo "APP: $app"
    if [ "$app" == "phmmer"  -o  "$app" == "jackhmmer" ]; then
	hmmer=1
	#appWrapper="$MDB_BENCHMARKING_SCRIPTS_DIR/hmmerWrapper.sh"
	pssmFile=""
    else
	hmmer=0
    fi
    appWrapper="$MDB_BENCHMARKING_SCRIPTS_DIR/${app}Wrapper.sh"
    exitIfFileDoesNotExist $appWrapper "ERROR: wrapper script \"$appWrapper\" does not exist (or is zero length)"

    for benchmarkDB in $DBS; do
	echo "DB: $benchmarkDB"

	if [ ! -d "$MDB_RUNS_BASE_DIR/$benchmarkDB" ]; then
	    cmd="mkdir $MDB_RUNS_BASE_DIR/$benchmarkDB"
	    if [ $VERBOSE -eq 1 ]; then echo $cmd; fi
	    $cmd
	fi
	submitDir="$MDB_RUNS_BASE_DIR/$benchmarkDB/$app$SUFFIX"
	if [ $DEBUG -eq 1 ]; then
	    echo "submitDir: $submitDir" 1>&2
	fi

        allDoneFile="$submitDir.allDone"
        if [ -f $allDoneFile ]; then
    	    if [ $VERBOSE -eq 1 ]; then echo "$submitDir is ALL DONE"; fi
    	    continue
        fi
	
	if [ ! -d $submitDir ]; then
	    cmd="mkdir $submitDir"
	    if [ $VERBOSE -eq 1 ]; then echo $cmd; fi
	    $cmd
	fi
	
	defaultQueryList=0
	if [ -z "$QUERIES" ]; then
	    if [ "$app" == "global" ]; then
		queriesFile="$MDB_DB_BASE_DIR/queriesList-msf.txt"
	    else
		queriesFile="$MDB_DB_BASE_DIR/$benchmarkDB/queriesList.txt"
	    fi
	    exitIfFileDoesNotExist $queriesFile
	    QUERIES=`cat $queriesFile`
	    if [ $DEBUG -eq 1 ]; then echo "QUERIES: $QUERIES"; fi
	    
	    defaultQueryList=1
	fi

	if [ "$app" == "global" ]; then
	    finalDatabase="$MDB_DB_BASE_DIR/finalDatabase-$app.fa"
	else
	    finalDatabase="$MDB_DB_BASE_DIR/finalDatabase.pin"
	fi
	relevanceFile="$MDB_DB_BASE_DIR/relevanceInfo.tab"
	for file in $finalDatabase $relevanceFile; do
	    exitIfFileDoesNotExist $file
	done

	allDone=1
			
	for queryFile in $QUERIES; do
	    if [ $DEBUG -eq 1 ]; then echo "queryFile: $queryFile"; fi
            queryFileBase=${queryFile##*/}
            queryFileBase=${queryFileBase%.*}
	    
            jobName="$queryFileBase.$app$SUFFIX"
            jobDir="$submitDir/$jobName"
            doneFile="$jobDir/.$jobName.done"
            
            if [ -f $doneFile ]; then
    		if [ $DEBUG -eq 1 ]; then echo "$jobName is DONE"; fi
    		continue
            fi
    	    
            inQueueTest=`perl -n -e "m/^$jobName\$/ && print" $jobsFile`
            if [ "$inQueueTest" ]; then
    	    	if [ $VERBOSE -eq 1 ]; then echo "$jobName is either pending or running"; fi
	    	allDone=0
    	    	continue
            fi
	    
	    if [ $ITERATIVE -eq 1 ]; then
		if [ $hmmer -eq 0 ]; then
		    pssmFile="$jobDir/${jobName}.nr.pssm"
		fi
	    fi
	    if [ $ITERATIVE -eq 0  -o  $PSSM_ONLY -ne 1 ]; then
		spougeFile="$jobDir/$jobName.spouge"
		if [ -n "$E_VALUE_THRESHOLD" ]; then
		    # spougeEFile initiated through $evalueStr
		    spougeEFile="$spougeFile.$E_VALUE_THRESHOLD"
		    spougeFile="$spougeEFile"
		fi
	    fi
            outFile="$jobDir/$jobName.out"
            errFile="$jobDir/$jobName.err"

	    #finalHitsFile="$jobDir/hits_$app.$queryFileBase.final.out"
	    finalHitsFile="$jobDir/$queryFileBase.$app.hits.final.txt"
	    outputFile="$jobDir/$queryFileBase.$app.final.out"
	    nrFile=""
	    if [ $ITERATIVE -eq 1 ]; then 
		#nrFile="$jobDir/$app.$queryFileBase.nr.out"
		#nrFile="$jobDir/$queryFileBase.$app.nr.out"
		nrFile="$jobDir/$queryFileBase.$app.iterationsDb.out"
	    fi
	    
            if [ ! -d $jobDir ]; then
    		if [ $VERBOSE -eq 1 ]; then
    		    echo "Directory \"$jobDir\" not found for $jobName";
    		fi
    		submitJob
    		continue
            fi
	    
            if [ ! -s $outFile  -a ! -f $outFile.gz ]; then
                echo "OUT file \"$outFile\" not found for $jobName";
    		submitJob
    		continue
            fi
	    
	    if [ $ITERATIVE -eq 1 ]; then 
		if [ ! -f $nrFile.gz ]; then
    		    echo "IterationDB file \"$nrFile.gz\" not found for $jobName";
    		    submitJob
    		    continue
		fi
            fi
	    
	    if [ $PSSM_ONLY -eq 1 ]; then
		if [ -n "$pssmFile"  -a  ! -s $pssmFile ]; then
    		    echo "PSSM file \"$pssmFile\" not found for $jobName";
    		    submitJob
    		    continue
		fi
	    else
		if [ -f "$errFile.gz" ]; then
		    unzippedErrFile=1
		    gunzip "$errFile.gz"
		else
		    unzippedErrFile=0
		fi
		
		if [ -s "$errFile" ]; then
		    # "No hits found in psisemiglobal.cd00007.iterationsDb.out! (skipping)"
    		    if [ `grep -l "No hits found" $errFile` ]; then
			grep "No hits found" $errFile | tee $doneFile
        		continue
		    fi
		    
    		    if [ `grep -l "ERROR: No successful subsets;" $errFile` ]; then
			grep "ERROR: No successful subsets;" $errFile | tee $doneFile
        		continue
		    fi
		fi	    
		
		spougeFileFound=1
		if [ ! -s $spougeFile ]; then
    		    echo "Spouge file \"$spougeFile\" not found for $jobName";
		    spougeFileFound=0
		fi

		if [ $spougeFileFound -eq 0 ]; then
   		    # test for "Sequence with id d1a8da2 no longer exists in database...alignment skipped"
		    if [ -s "$finalHitsFile" ]; then
			grepResult=`grep "no longer exists in database" $finalHitsFile`
			if [ -n "$grepResult" ]; then
			    echo "    Final hits file \"$finalHitsFile\" has \"no longer exists in database\"!  Removing."
			    rm $finalHitsFile
			fi
		    fi
    	     	    submitJob
    	     	    continue
		fi
    		
		if [ $ROCN_SPOUGE_FILES -eq 1  -a   ! -s "${spougeFile%.spouge}-rocn.spouge" ]; then
    		    echo "ROC_n Spouge file \"${spougeFile%.spouge}-rocn.spouge\" not found for $jobName";
    		    submitJob
    		    continue
		fi
		
		if [ -f "$outFile.gz" ]; then
		    unzippedOutFile=1
		    gunzip "$outFile.gz"
		else
		    unzippedOutFile=0
		fi
		
		if [ ! `grep -l "$FINISHED_LINE" $outFile` ]; then
    		    echo "OUT file \"$outFile\" not complete for $jobName";
    		    submitJob
    		    continue
		fi


		if [ $unzippedOutFile -eq 1 ]; then
		    gzip "$outFile"
		fi
		
		if [ $unzippedErrFile -eq 1 ]; then
		    gzip "$errFile"
		fi
	    fi # END PSSM_ONLY
	    
	    
        #if [ $VERBOSE -eq 1 ]; then
            echo "$jobName is now DONE";
        #fi
            
            touch $doneFile

	done # END QUERIES

	if [ $AGGREGATE_JOBS -eq 1  -a  $AGGREGATE_JOBS_COUNT -gt 0  -a $AGGREGATE_JOBS_COUNT -ne $AGGREGATE_JOBS_NUMBER ]; then
	    submitAggregateJobs
	fi


	if [ $allDone -eq 1  -a  $defaultQueryList -eq 1 ]; then
            echo "$submitDir is now ALL DONE";
	    touch $allDoneFile
	fi
    done     # END DBS
done         # END APPS

cleanUpAndExit 0
