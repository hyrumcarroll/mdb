#!/bin/bash
	
#$ -N TIMING_TEMPLATE
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -j n
## -o /panfs/pan1/iglobal/projs/multiDomainBenchmarkRuns/MultiDomainBenchmark-test-timing/TIMING_TEMPLATE.out.job
## -e /panfs/pan1/iglobal/projs/multiDomainBenchmarkRuns/MultiDomainBenchmark-test-timing/TIMING_TEMPLATE.err.job
#$ -P unified
#$ -l h_rt=360000
#$ -l h_vmem=8.1G
#$ -l mem_free=6G
# use reserve_mem only for large jobs b/c it will block other users from getting access to that node
#$ -l reserve_mem=8G
## -m a
#$ -m n


source ~/.bashrc
ulimit -v 8388608

hostname
cat /proc/cpuinfo
cat /proc/meminfo

names="NAMES_TEMPLATE"
apps="jackhmmer phmmer psiblast psiblast-non_iterative psisemiglobal psisemiglobal-non_iterative"
for name in $names; do
    echo "name: $name"
    for app in $apps; do
	echo "app: $app"
	appDir="/panfs/pan1/iglobal/projs/multiDomainBenchmarkRuns/MultiDomainBenchmark-test-timing/TIMING_TEMPLATE/$app"
	
	jobName="$name.$app"
	queryFile="/panfs/pan1.be-md.ncbi.nlm.nih.gov/iglobal/dbs/multiDomainBenchmark/queries-multiDomain/$name.fa"
	spougeFile="$appDir/$name.$app/$name.$app.spouge.1000"
	outFile="$appDir/$name.$app/$name.$app.out"
	errFile="$appDir/$name.$app/$name.$app.err"
	iterative=1
	hitsFile="$appDir/$name.$app/$name.$app.hits.final.txt"
	outputFile="$appDir/$name.$app/$name.$app.final.out"
	if [ $iterative -eq 1 ]; then 
	    pssmFile="$appDir/$name.$app/$name.$app.pssm"
	    if [ -n "$pssmFile" ]; then
		pssmStr=" -pssm=$pssmFile "
		pssmOnly=0
		if [ $pssmOnly -eq 1 ]; then
		    pssmStr=" $pssmStr -pssmOnly "
		fi
	    fi
	    iterativeStr=" $pssmStr   "
	fi

	if [ ! -d "$appDir" ]; then 
	    mkdir -p "$appDir"
	fi
	cd "$appDir"

	if [ ! -d "$jobName" ]; then
	    cmd="mkdir $jobName"
    #echo "$cmd"
	    $cmd
	fi
	cd $jobName

	for filename in $outFile $errFile; do
	    if [ -f $filename.gz ]; then
		gunzip $filename.gz
	    fi
	    echo "-----------------------------------------------------------------------------"`date` >> $filename
	done    

	flagsFileStr="";
	flagsFilename="$HOME/research/mdb/benchmarkingScripts/flagsFile_${app}_.txt"
	
	if [ -s "$flagsFilename" ]; then
	    #flagsFileStr="-flagsFile="`cat "$flagsFilename"`
	    flagsFileStr="-flagsFile=$flagsFilename"
	fi
	
	flagsFinalFileStr="";
	flagsFinalFilename="$HOME/research/mdb/benchmarkingScripts/flagsFinalFile_${app}_.txt"

	if [ -s "$flagsFinalFilename" ]; then
	    #flagsFinalFileStr="-flagsFinalFile="`cat "$flagsFinalFilename"`
	    flagsFinalFileStr="-flagsFinalFile=$flagsFinalFilename"
	fi
	
	cmd="$HOME/research/MultiDomainBenchmark/benchmarkingScripts/${app}Wrapper.sh -$app -query=$queryFile  -spouge=$spougeFile -db=/panfs/pan1/iglobal/dbs/multiDomainBenchmark/finalDatabase.pin  -rel=/panfs/pan1/iglobal/dbs/multiDomainBenchmark/relevanceInfo.tab   -v  $flagsStr  $flagsFinalFileStr  -hits=$hitsFile -out=$outputFile   -evalue=1000   $iterativeStr "

	echo "$cmd" >> $outFile 2>> $errFile
	if [ 0 -eq 1 ]; then exit; fi
	$cmd >> $outFile 2>> $errFile

	rc=$?
	if [ $rc != 0 ]; then
	    echo ":$cmd: FAILED (rc = $rc)" >> $errFile
	    exit 1
	fi

	if [ -n "$pssmFile" ]; then
	    if [ $pssmOnly -eq 1 ]; then
		if [ ! -s $pssmFile ]; then
		    echo "PSSM file, $pssmFile, does not exist!" >> $errFile
		    exit 1
		fi
	    elif [ ! -s $pssmFile.gz ]; then
		echo "PSSM file, $pssmFile.gz, does not exist!" >> $errFile
		exit 1
	    fi
	fi

	if [ ! -s $spougeFile ]; then
	    echo "Spouge file, $spougeFile, does not exist!" >> $errFile
	    exit 1
	fi

	if [ ! -s $errFile ]; then
	    rm $errFile
	fi

	cd "$appDir"

	(echo ""; echo "Finished - ALL DONE") >> $outFile
	echo ""

    done # apps
    echo ""
done # names

exit 0
