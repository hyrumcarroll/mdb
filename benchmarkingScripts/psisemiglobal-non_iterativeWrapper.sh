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

sh $MDB_BENCHMARKING_SCRIPTS_DIR/psi-non_iterativeWrapper.sh  -psisemiglobal-non_iterative  $@

exit 0
