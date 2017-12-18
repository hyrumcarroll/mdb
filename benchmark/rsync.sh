#!/bin/bash

backupFile.pl ~/research/mdb/benchmark/*.sh ~/research/mdb/benchmark/*.pl ~/research/mdb/benchmark/bak/

ssh $ranger "rsync -av \$air:research/mdb/benchmark/ research/mdb/benchmark/" && rsync -av $ranger:~/research/mdb/benchmark/ ~/research/mdb/benchmark/
