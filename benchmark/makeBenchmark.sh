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

# Making multi-domain queries:
# 1. Simplify the domain locations by ignoring domains that are embedded within themselves
# 1.a. Remove the location annotation and remove the sequence from the library for the embedded domains
# 2. Compile domain architectures (DAs) for each library sequence
# 3. Enumerate entries that have multiple domains
# 4. Choose random multiple domain sequences as queries

# If REFPROTDOM_DIR is not an environment variable, then default to the current directory
if [ -z "$REFPROTDOM_DIR" ]; then
    REFPROTDOM_DIR="."
fi

VERBOSE=1
DEBUG=0

CALCULATE_STATISTICS=1
USE_CACHED_FILES=1
JUST_REMOVE_PRODUCTS=0

verbosityStr=""

function VERBOSE {
    if [ $VERBOSE -ge 1 ]; then
	echo "$@"
    fi
}

function DEBUG {
    if [ $DEBUG -ge 1 ]; then
	echo "DEBUG: $@"
    fi
}

function run {
    cmd=$1
    shift 1
    products="$@"

    if [ $DEBUG -ge 2 ]; then 
	DEBUG "run()"
	DEBUG "    cmd=\"$cmd\""
	DEBUG "    products=\"$products\""
	time -p eval "$cmd"
    else
	VERBOSE "$cmd"
	eval "$cmd"
    fi
    
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


usageStr="Usage $0:
  [-cache]            Use files produced from earlier executions (default: set)
  [-nocache]          Do not use files produced from earlier executions (default: unset)
  [-stats]            Calculate summary statistics (default: unset)
  [-clean]            Only remove files produced by this script (default: unset)
  [-refprotdom dir]   Directory of RefProtDom v1.2 files (namely, family_members.annot
                      and library_all_domains.fa) (defaults to environment variable
                      REFPROTDOM, if it's set, otherwise the current directory)
  [-v | -verbose]     Echo commands and other information (additonal -d arguments
                      increases the amount of information displayed) (default: set)
  [-nov | -noverbose] Silence most of the output (default: unset)
  [-d | -debug]       Display debugging information (additonal -d arguments increases
                      the amount of information displayed) (default: unset)
"         

USE_CACHED_FILES_ARG=1
HELP=0

while [ $# -ge 1 ]; do
    arg="$1"
    shift

    case "$arg" in
	-v|-verbose|--verbose         )  VERBOSE=$(( VERBOSE + 1 )) ;;
	-nov|-noverbose|--noverbose   )  VERBOSE=0 ;;
	-d|-debug|--debug             )  DEBUG=$(( DEBUG + 1 )) ;; 

	-cache|--cache                )  USE_CACHED_FILES=1 ;;
	-nocache|--nocache            )  USE_CACHED_FILES=0 ;;
	
	-stat|-stats|--stat|--stats   )  CALCULATE_STATISTICS=1 ;;

	-refprotdom|--refprotdom      )  REFPROTDOM_DIR="$1"; shift ;;

	-h|-help|--help               )  HELP=1 ;;
	-clean|--clean                )  JUST_REMOVE_PRODUCTS=1 ;;
		       
	*                             )  echo "ERROR: Unknown arg: $arg!" >&2; exit 1;;
    esac
done 

verbosityStr=""
# pass along -v and -d, but decremented by one level
i=2
while [ $i -le $VERBOSE ]; do
    verbosityStr="$verbosityStr -v"
    echo "verbosityStr (-v): $verbosityStr"
    i=$(( i + 1 ))
done

i=2
while [ $i -le $DEBUG ]; do
    verbosityStr="$verbosityStr -d"
    echo "verbosityStr (-d): $verbosityStr"
    i=$(( i + 1 ))
done

if [ $DEBUG -ge 1 ]; then
    echo "$0 parameters:"
    echo "USE_CACHED_FILES:      $USE_CACHED_FILES"
    echo "CALCULATE_STATISTICS:  $CALCULATE_STATISTICS"
    echo "JUST_REMOVE_PRODUCTS:  $JUST_REMOVE_PRODUCTS"
    echo "REFPROTDOM_DIR:        $REFPROTDOM_DIR"
    echo "VERBOSE:               $VERBOSE"
    echo "DEBUG:                 $DEBUG"
    echo "verbosityStr:          $verbosityStr"
    echo ""
fi

if [ $HELP -eq 1 ]; then
    echo "$usageStr"
    exit 0
fi

scriptsDir="scripts"

nonOverlappingDomainLocsFilename="family_members-nonOverlapping.annot"
multiDomainFullLabelsOutputFilename="domainQueries_multi.tab"
singleDomainFullLabelsOutputFilename="domainQueries_single.txt"

taxaRemoved="sequencesWithOverlappingDomains.txt"
libraryFilename="library_all_domains_MDB.fa"
libraryTaxonAndSpeciesCodesFilename="${libraryFilename%.fa}-taxa_speciesCodes.tab"

# in the taxonomy directory
taxonomyFilename="taxonomy/mdbTaxonomy.tab"  
speciesListFilename="taxonomy/speclist.txt" # cache file
olderSpeciesListFilename="taxonomy/speclist-20060516.txt"
daDiversityFilename="taxonomy/dasDiversity.tab"
diverseDAsFilename="taxonomy/diverseDAs.tab"
taxonomyErrorsFilename="taxonomy/taxonomy.err"
uniprotMappingsCachedFile="taxonomy/idmapping_selected-cached.tab" # cache file
ncbiTaxonomyNamesCacheFile="taxonomy/names.dmp"
ncbiTaxonomyNodesCacheFile="taxonomy/nodes.dmp"
ncbiTaxonomyBaseFilename="taxonomy/taxdump" #.tar.gz

getDAsErrorFilename="getDAs.pl.err"
getDAsOutFilename="getDAs.pl.out"
daInfoFilename="daInfo.tab"
relevanceInfoFilename="taxon2da.tab"
allQueriesFilename="queriesPerDA-all.tab"
diverseQueriesFilename="queriesPerDA-diverse.tab"

domainLocsFilename="$nonOverlappingDomainLocsFilename"
multiDomainOutputFilename="$multiDomainFullLabelsOutputFilename"
singleDomainOutputFilename="$singleDomainFullLabelsOutputFilename"


domainLocsSimplifiedFilename="family_members-nonOverlapping-simplified.annot"

dasPerQueryHistogramFilename="numDAsPerQuery.histo"
daSizesFilename="daSizes.csv"
numDomainsHistogramFilename="numDomains-queries.histo"
domainSizesFilename="$multiDomainOutputFilename.queries"

diverseQueriesTrainingFilename="${diverseQueriesFilename%.tab}-training.tab"
diverseQueriesTestFilename="${diverseQueriesFilename%.tab}-test.tab"

fastaBaseName="${libraryFilename%.fa}"
classifyCriterionFilename="$fastaBaseName.classifyCriterion"

if [ $USE_CACHED_FILES == 0  -o  $JUST_REMOVE_PRODUCTS == 1 ]; then 
    # remove all products
    ALL_PRODUCTS="$relevanceInfoFilename  $multiDomainFullLabelsOutputFilename $domainSizesFilename $singleDomainFullLabelsOutputFilename   $nonOverlappingDomainLocsFilename  formatdb.log  $daSizesFilename   lengths-test.tab  lengths-training.tab  lengths.tab  $libraryFilename  longestQueriesForDAs.txt  $dasPerQueryHistogramFilename  $numDomainsHistogramFilename  queries-multiDomain/  queriesList-test.txt  queriesList-training.txt   $allQueriesFilename  $diverseQueriesFilename $diverseQueriesTrainingFilename $diverseQueriesTestFilename   $taxaRemoved  shortestQueriesForDAs.txt $domainLocsSimplifiedFilename  MultiDomainBenchmark-training  MultiDomainBenchmark-test  $fastaBaseName.{classifyCriterion,phr,pin,psq}  $classifyCriterionFilename  finalDatabase.{classifyCriterion,fa,phr,pin,psq} domainLocs.tab  $libraryTaxonAndSpeciesCodesFilename $speciesListFilename $daInfoFilename $daDiversityFilename $diverseDAsFilename $getDAsErrorFilename $getDAsOutFilename $taxonomyFilename $ncbiTaxonomyNamesCacheFile{,.gz} $ncbiTaxonomyNodesCacheFile{,.gz}  $ncbiTaxonomyBaseFilename{,.tar,.tar.gz} $taxonomyErrorsFilename"
    # makeBenchmark.sh.out
    run "rm -fr $ALL_PRODUCTS 2> /dev/null"
    #run "rm -f $uniprotMappingsCachedFile 2> /dev/null"

    if [ $JUST_REMOVE_PRODUCTS == 1 ]; then
	exit 0
    fi
fi


scripts="$scriptsDir/domainLocsParser.pl  $scriptsDir/getDAs.pl  $scriptsDir/removeTaxa.pl  $scriptsDir/fileSplicer.pl"
if [ $CALCULATE_STATISTICS -eq 1 ]; then
    scripts="$scripts  $scriptsDir/getDASizes.pl  $scriptsDir/getDomainSizes.pl  $scriptsDir/lengths.sh"
fi
for script in $scripts; do
    if [ ! -e "$script" ]; then
	echo "ERROR: Script \"$script\" is not executable!" >&2
	exit 1
    fi
done

# 1. Simplify the domain locations by ignoring domains that are embedded within themselves
if [ ! -s "$nonOverlappingDomainLocsFilename"  -o  ! -s "$multiDomainFullLabelsOutputFilename"  -o  ! -s "$singleDomainFullLabelsOutputFilename"  -o ! -s "$taxaRemoved"  -o  ! -s "$libraryFilename" -o ! -s "$taxonomyFilename" ]; then
    domainLocsOriginalFilename="$REFPROTDOM_DIR/family_members.annot"
    if [ ! -s "$domainLocsOriginalFilename" ]; then
	echo "ERROR: Required RefProtDom input file \"$domainLocsOriginalFilename\" not found! (Note: REFPROTDOM_DIR is set to \"$REFPROTDOM_DIR\")" >&2
	exit 1
    fi

    # create (temporary) copy of RefProtDom domain locs file and replace '|'s with '_'s
    domainLocsTempFilename="family_members.annot.tmp"
    if [ ! -s "$domainLocsTempFilename" ]; then
	cmd="perl -p -e 's/\|/_/g'  < $domainLocsOriginalFilename  > $domainLocsTempFilename"
	run "$cmd" "$domainLocsTempFilename"
    fi
    
    cmd="$scriptsDir/domainLocsParser.pl $verbosityStr -in $domainLocsTempFilename -out $nonOverlappingDomainLocsFilename -multi $multiDomainFullLabelsOutputFilename -single $singleDomainFullLabelsOutputFilename -removed $taxaRemoved"
    run "$cmd"  $nonOverlappingDomainLocsFilename $multiDomainFullLabelsOutputFilename $singleDomainFullLabelsOutputFilename $taxaRemoved

    # remove temp file
    run "rm $domainLocsTempFilename"

    
    libraryOriginalFilename="$REFPROTDOM_DIR/library_all_domains.fa"
    if [ ! -s "$libraryOriginalFilename" ]; then
	echo "ERROR: Required RefProtDom input file \"$libraryOriginalFilename\" not found! (Note: REFPROTDOM_DIR is set to \"$REFPROTDOM_DIR\")" >&2
	exit 1
    fi

    # create (temporary) copy of RefProtDom FASTA file and replace '|'s with '_'s
    libraryTempFilename="library_all_domains.fa.tmp"
    if [ ! -s "$libraryTempFilename" ]; then
	cmd="perl -p -e 's/\|/_/g'  < $libraryOriginalFilename  > $libraryTempFilename"
	run "$cmd" "$libraryTempFilename"
    fi

    # run "wc -l $libraryTempFilename $taxaRemoved"  # DEBUGGING
    cmd="$scriptsDir/removeTaxa.pl $verbosityStr -in $libraryTempFilename -out $libraryFilename -removeFile $taxaRemoved"
    run "$cmd"  "$libraryFilename"
    #run "wc -l $libraryTempFilename $libraryFilename $taxaRemoved"  # DEBUGGING

    # remove temp file
    run "rm $libraryTempFilename"

    #
    # Parse the MDB .fa for taxon labels and split them up into <accession>\t<species_code> tuples
    #
    cmd="grep '^>' $libraryFilename > $libraryTaxonAndSpeciesCodesFilename"
    run "$cmd"  "$libraryTaxonAndSpeciesCodesFilename"

    # e.g.: >pfam21_O00098_CISY_EMENI
    # Parse out:    ^^^^^^ and  ^^^^^
    perl -pi -e 's/^>[^_]+_([^_]+).*_([^_]+)/$1\t$2/' "$libraryTaxonAndSpeciesCodesFilename"
    
    # run "head $libraryTaxonAndSpeciesCodesFilename"

    #
    # Get taxonomy information to inform the query selection
    #
    # cmd="cd taxonomy; ./uniprot2taxonomy.pl  --uniprots $libraryTaxonAndSpeciesCodesFilename  --out $taxonomyFilename  --speclist $speciesListFilename  --speclist2 $olderSpeciesListFilename  --uniprotMappings $uniprotMappingsCachedFile  2> $taxonomyErrorsFilename; cd ../"
    cmd="./taxonomy/uniprot2taxonomy.pl  $verbosityStr  --uniprots $libraryTaxonAndSpeciesCodesFilename  --out $taxonomyFilename  --speclist $speciesListFilename  --speclist2 $olderSpeciesListFilename  --names $ncbiTaxonomyNamesCacheFile  --nodes $ncbiTaxonomyNodesCacheFile  --uniprotMappings $uniprotMappingsCachedFile --tempDir taxonomy  2> $taxonomyErrorsFilename"
    run "$cmd"  "$taxonomyFilename"

    run  "grep -v '^ALERT:' $taxonomyErrorsFilename"
else
    VERBOSE "$nonOverlappingDomainLocsFilename, $multiDomainFullLabelsOutputFilename, $singleDomainFullLabelsOutputFilename, $taxaRemoved, $libraryFilename, and $taxonomyFilename already exist"
fi


if [ ! -s "$relevanceInfoFilename"  -o ! -s "$allQueriesFilename"  -o ! -s "$daInfoFilename"  -o ! -s "$diverseDAsFilename" ]; then
    # getDAs.pl
    #   Determines all of the domain architectures (DAs) and which sequences belong to each one.
    #   Also filters potential queries with global variables: $MIN_QUERY_LENGTH = 10; $MAX_QUERY_LENGTH = 1800;
    products="$relevanceInfoFilename  $allQueriesFilename  $daInfoFilename"
    cmd="$scriptsDir/getDAs.pl  $verbosityStr  --domainLocs $nonOverlappingDomainLocsFilename  --fasta $libraryFilename  --rel $relevanceInfoFilename  --queries $allQueriesFilename  --da $daInfoFilename  > $getDAsOutFilename  2> $getDAsErrorFilename"
    run "$cmd" $products

    #
    # Determine which domains the members of each of the DAs belong to (for those DAs that are candidates for queries)
    #
    #cmd="cd taxonomy; ./diverseDAs.py  ../$daInfoFilename  $taxonomyFilename ../$allQueriesFilename > $daDiversityFilename; cd ../"
    cmd="./taxonomy/diverseDAs.py  $daInfoFilename  $taxonomyFilename $allQueriesFilename  $verbosityStr > $daDiversityFilename"
    run "$cmd"  "$daDiversityFilename"

    #
    # Get DAs that have members in 2 or 3 domains (that are candidates for queries)
    #
    cmd="head -n 1 $daDiversityFilename > $diverseDAsFilename && perl -n -e 'm/^(\S+)\s[23]/ && print' $daDiversityFilename >> $diverseDAsFilename"
    run "$cmd" "$diverseDAsFilename"
    
    # run "wc -l $diverseDAsFilename"
else
    VERBOSE "$relevanceInfoFilename, $allQueriesFilename, $daInfoFilename, and $diverseDAsFilename already exist"
fi


if [ ! -s "$diverseDAsFilename"  -o  "$allQueriesFilename" -ot "$diverseDAsFilename"  -o "$diverseDAsFilename" -ot "$diverseDAsFilename" ]; then
    #
    # Randomly select diverse query sequences
    #
    # get a list of queries in $allQueriesFilename filtered by $diverseDAsFilename
    cmd="$scriptsDir/filterByUnion.pl  $verbosityStr -main $allQueriesFilename  -keys $diverseDAsFilename  >  $diverseQueriesFilename"
    run "$cmd" "$diverseQueriesFilename"
fi

if [ ! -s "$domainLocsSimplifiedFilename" ]; then
    #echo "Making $domainLocsFilename"
    cmd="perl -p -e 's/^(\S+\s+\S+\s+\S+).+/\$1/' $domainLocsFilename > $domainLocsSimplifiedFilename"
    run "$cmd" "$domainLocsSimplifiedFilename"
else
    VERBOSE "$domainLocsSimplifiedFilename already exists"
fi

if [ $CALCULATE_STATISTICS -eq 1 ]; then
    if [ ! -s "$dasPerQueryHistogramFilename"  -o  "$dasPerQueryHistogramFilename" -ot "$diverseQueriesFilename"  -o  "$dasPerQueryHistogramFilename" -ot "$relevanceInfoFilename" ]; then
	if [ ! -s "$daSizesFilename"  -o  "$daSizesFilename" -ot "$diverseQueriesFilename"  -o  "$daSizesFilename" -ot "$relevanceInfoFilename" ]; then
	    # prints out "taxon,DA,DA_size"
	    cmd="$scriptsDir/getDASizes.pl $diverseQueriesFilename $relevanceInfoFilename > \"$daSizesFilename\""
	    run "$cmd" "$daSizesFilename"
	fi

	echo "Making histogram of the number of DAs for each query . . ."
	echo "# Number_of_DAs DA_Membership_count" >  "$dasPerQueryHistogramFilename"
	cmd="grep -v \"^#\" $daSizesFilename | cut -f 2 -d ',' | sort -n | uniq -c  >>  \"$dasPerQueryHistogramFilename\""
	run "$cmd"  "$dasPerQueryHistogramFilename"
    fi


    if [ ! -s "$numDomainsHistogramFilename"  -o  "$numDomainsHistogramFilename" -ot "$diverseQueriesFilename"  -o  "$numDomainsHistogramFilename" -ot "$multiDomainOutputFilename" ]; then
	if [ ! -s "$domainSizesFilename"  -o  "$domainSizesFilename" -ot "$diverseQueriesFilename"  -o  "$domainSizesFilename" -ot "$multiDomainOutputFilename" ]; then
	    cmd="$scriptsDir/getDomainSizes.pl $diverseQueriesFilename  $multiDomainOutputFilename   >  \"$domainSizesFilename\""
	    run "$cmd"  "$domainSizesFilename"
	fi

	echo "Making histogram of the number of domains in the random queries . . ."
	echo "# Number_of_DAs Number_of_Domains" >  "$numDomainsHistogramFilename"
	cmd="grep -v \"^#\" \"$domainSizesFilename\" | cut -f 2 | sort -n | uniq -c  >>  \"$numDomainsHistogramFilename\""
	run "$cmd" "$numDomainsHistogramFilename"
    fi
fi

if [ ! -s "$diverseQueriesTrainingFilename"  -o  ! -s "$diverseQueriesTestFilename"  -o  "$diverseQueriesTrainingFilename" -ot "$diverseQueriesFilename"  -o  "$diverseQueriesTestFilename" -ot "$diverseQueriesFilename" ]; then
    
    echo "Spliting queries into training and test sets . . ."
    grep -v "^#" $diverseQueriesFilename > $diverseQueriesFilename.tmp # remove header line
    cmd="$scriptsDir/fileSplicer.pl $diverseQueriesFilename.tmp 2"
    run "$cmd"  $diverseQueriesFilename.tmp.0 $diverseQueriesFilename.tmp.1
    rm $diverseQueriesFilename.tmp

    mv -f  $diverseQueriesFilename.tmp.0   "$diverseQueriesTrainingFilename"
    mv -f  $diverseQueriesFilename.tmp.1   "$diverseQueriesTestFilename"

    
    queriesDir="queries-multiDomain"
    if [ ! -d "$queriesDir" ]; then
	mkdir "$queriesDir"
    fi
    for testType in training test; do 
	echo "Making/getting multi-domain $testType queries . . ."
	queriesListFilename="queriesList-$testType.txt"
	if [ -s "$queriesListFilename" ]; then
	    if [ -e "backupFile.pl" ]; then
		backupFile.pl "$queriesListFilename"
	    fi
	    rm -f "$queriesListFilename"
	fi
	# # replaced the following loop with a perl script for efficency 
	# for taxon in `cat ${diverseQueriesFilename%.tab}-$testType.tab`; do
	#     taxonFilename="${taxon//|/_}"  # replace all occurrences of "|" with "_"
	#     fastaFilename="$queriesDir/$taxonFilename.fa"
	#     grep -A 1 "^>$taxon"  "$libraryFilename" > "$fastaFilename"
	#     echo "$fastaFilename" >> $queriesListFilename
	# done

	# make labels file, then extend it to be filenames for each of the queries
	cmd="cut -f 2 ${diverseQueriesFilename%.tab}-$testType.tab > $queriesListFilename"
	run "$cmd"  "$queriesListFilename"
	
	run "$scriptsDir/splitDataset.pl -in $libraryFilename  -labels $queriesListFilename   -dir $queriesDir"
	
	# echo "Compiling $queriesListFilename . . ."
	# for taxon in `cat ${diverseQueriesFilename%.tab}-$testType.tab`; do
	#     taxonFilename="${taxon//|/_}"  # replace all occurrences of "|" with "_"
	#     fastaFilename="$queriesDir/$taxonFilename.fa"
	#     echo "$fastaFilename" >> $queriesListFilename
	# done

	# changes labels (in $queriesListFilename) to be full path filenames to each query
	perl -pi -e "s%\|%_%g" $queriesListFilename  # replace all occurrences of "|" with "_"
	perl -pi -e "s%^%$PWD/$queriesDir/%" $queriesListFilename  # make the query filenames full paths:
	perl -pi -e 's%$%.fa%' $queriesListFilename  # make the query filenames full paths:
	
	#cmd="wc -l $queriesListFilename"
	#run "$cmd"
    done
else
    VERBOSE "$queriesTrainingFilename and $queriesTestFilename already exist"
fi

if [ ! -L "domainLocs.tab" ]; then 
    cmd="ln -s $domainLocsSimplifiedFilename  domainLocs.tab"
    run "$cmd" "domainLocs.tab"
fi


# make BLAST DBs
if [ `which "makeblastdb"` ]; then
    echo "Making BLAST DBs . . ."
    blastDbFilename="${libraryFilename%.fa}.pin"
    if [ ! -s "$blastDbFilename"  -o  "$libraryFilename" -nt "$blastDbFilename" ]; then
	# make BLAST DBs
	cmd="makeblastdb -in $libraryFilename -dbtype prot -out $fastaBaseName"
	run "$cmd" $fastaBaseName.{phr,pin,psq}
	#ls -l library_all_domains.p*
    fi
    if [ ! -s "$classifyCriterionFilename" ]; then
	echo " --domainLocs=$PWD/domainLocs.tab  --overlap=50  --norandomsAsIrrelevants  --combineHSPs" > "$classifyCriterionFilename"
    fi

    # make links for pipeline
    for i in classifyCriterion fa phr pin psq; do
	if [ ! -L "finalDatabase.$i" ]; then 
	    cmd="ln -s ${libraryFilename%.fa}.$i finalDatabase.$i"
	    run "$cmd" "finalDatabase.$i"
	fi
    done
else
    echo "NOTE: makeblastdb not found (skipping making the BLAST databases)" 
fi

if [ ! -L "relevanceInfo.tab" ]; then 
    cmd="ln -s $relevanceInfoFilename  relevanceInfo.tab"
    run "$cmd" "relevanceInfo.tab"
fi


for testType in training test; do 
    echo "Making $testType links . . ."
    queriesListFilename="queriesList-$testType.txt"

    dirName="MultiDomainBenchmark-$testType"
    if [ ! -d "$dirName" ]; then
	run "mkdir $dirName" "$dirName"
    fi
    
    cd "$dirName"
    queriesListFilenameLocal="queriesList.txt"
    if [ ! -L "$queriesListFilenameLocal" ]; then
	run "ln -s ../$queriesListFilename  $queriesListFilenameLocal" "$queriesListFilenameLocal"
    fi
    if [ ! -L "queriesDir" ]; then
	run "ln -s ../queries-multiDomain  queriesDir" "queriesDir"
    fi
    
    for filename in domainLocs.tab  finalDatabase.{classifyCriterion,fa,phr,pin,psq}  $relevanceInfoFilename  relevanceInfo.tab; do
	linkFilename="../$filename"
	if [ ! -L "$filename" ]; then 
	    if [ ! -s "../$filename" ]; then
		echo "ALERT: Creating a link in $dirName to ../$filename EVEN THOUGH IT DOES NOT EXIST!" >&2
	    fi
	    run "ln -s ../$filename  ."
	fi
    done
    
    cd -
done
	

if [ $CALCULATE_STATISTICS -eq 1 ]; then
    run "$scriptsDir/lengths.sh"
fi

echo "All done!"

exit 0
