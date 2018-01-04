#!/bin/bash

#
# randomly selected queries
#
names=`cat $MDB_RUNS_BASE_DIR/MultiDomainBenchmark-test/sampledDatasets.txt`
numReplicates=11

templateScriptFilename="timing-template.sh"

for i in `seq -w 1 $numReplicates`; do
    timingScriptFilename="timing$i.sh"
    perl -p  -e "s/TIMING_TEMPLATE/timing$i/g" "$templateScriptFilename" > "$timingScriptFilename"
    perl -pi -e "s/NAMES_TEMPLATE/$names/" "$timingScriptFilename"
    chmod 755 "$timingScriptFilename"
    echo qsub "$timingScriptFilename"
done

qstat
