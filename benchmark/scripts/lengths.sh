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

for subset in "training" "test"; do
    outputFileName="lengths-$subset.tab"
    if [ -s "$outputFileName"  -a -e "backupFile.pl" ]; then
	backupFile.pl "$outputFileName" 
    fi

    echo "#number_residues	query_file" > "$outputFileName"

    queriesList="queriesList-$subset.txt"
    for queryFile in `cat $queriesList`; do
	if [ ! -s "$queryFile" ]; then
	    echo "ERROR: queryFile: $queryFile" >&2
	    exit 1
	fi
	#echo "queryFile: $queryFile"
	
	numResidues=`head -n 2 "$queryFile" | tail -n 1 | wc -c`
	numResidues=$((numResidues - 1))  # to account for the newline
	echo "$numResidues	$queryFile"
    done | sort -n  >> "$outputFileName"
done

# combine subset files into single file

outputFileName="lengths.tab"
if [ -s "$outputFileName"  -a -e "backupFile.pl" ]; then
    backupFile.pl "$outputFileName" 
fi

echo "#number_residues	query_file" > "$outputFileName"
grep -v "^#" lengths-*.tab | sort -n >> "$outputFileName"
