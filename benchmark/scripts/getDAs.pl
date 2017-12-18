#!/usr/bin/env perl

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

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

BEGIN
{
    # Look for perl modules in ./perl
    unshift( @INC, "./perl"); 
}

use HDC::Common qw(:ALL);
use HDC::SortingShuffling qw( &fisherYatesShuffle);

my $MIN_QUERY_LENGTH = 10;
my $MAX_QUERY_LENGTH = 1800;

my $NUMBER_OF_DOMAIN_FIELDS = 6;  # Number of columns in the domain locations input file (from RefProtDom)
my $DOMAIN_NAME_FIELD = 0;        # Assumes domain names are the first field
my $MODE_FIELD = 0;               # Assumes domain names are the fourth field

# input

my $fastaFilename = "library_all_domains.fa";
my $domainLocsFilename = "family_members-nonOverlapping.annot";

# output
my $daInformationFilename = "daInfo.tab";
my $relevanceFilename = "taxon2da.tab";
my $queriesFilename   = "queriesPerDA.tab";

my $man = 0;
my $help = 0;

unless( GetOptions("-v|verbose+"        => \$VERBOSE,
		   "-d|debug+"          => \$DEBUG,
		   "fa|fas|fasta=s"     => \$fastaFilename,
		   "domainLocs=s"       => \$domainLocsFilename,
		   "da=s"               => \$daInformationFilename,
		   "rel=s"              => \$relevanceFilename,
		   "queries=s"          => \$queriesFilename,
		   'help|?'             => \$help,
		   'man'                => \$man)){
    pod2usage(2);    
}

DEBUG( "$0 parameters:
    fastaFilename:                $fastaFilename
    domainLocsFilename:           $domainLocsFilename
    daInformationFilename:        $daInformationFilename
    relevanceFilename:            $relevanceFilename
    queriesFilename:              $queriesFilename
    help:                         $help
    man:                          $man
    VERBOSE:                      $VERBOSE
    DEBUG:                        $DEBUG
");

pod2usage(-exitstatus => 0, -verbose => 1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

if( ! -s $fastaFilename){
    die( "ERROR: Fasta (library) filename \"$fastaFilename\" does not exist: $!");
}

if( ! -s $domainLocsFilename){
    die( "ERROR: Family filename \"$domainLocsFilename\" does not exist: $!");
}

my $fastaFileContents = fileToString( $fastaFilename);
my @seqs = split( '>', $fastaFileContents);
my %seqLens; # keys: taxon names / GI; values: sequence lengths
print STDERR "DEBUG: Found " . ($#seqs + 1 - 1) ." seqs in $fastaFilename\n";
if( $seqs[0] !~ m/^\s*$/){
    die( "ERROR: seqs[0]: $seqs[0]");
}
# ignore the first result from split (so seqIndex starts at 1)
for( my $seqIndex = 1; $seqIndex <= $#seqs; $seqIndex++){
    my ($taxon, $seq) = split( /\n/, $seqs[$seqIndex]);
    if( ! defined( $taxon)  || ! defined( $seq)){
	die( "ERROR: Could not parse index $seqIndex from $fastaFilename: $seqs[$seqIndex]\n");
    }
    $seqLens{$taxon} = length($seq);
}

my $domainLocsFileContents = fileToString( $domainLocsFilename);
#print "domainLocsFileContents:\n" . substr( $domainLocsFileContents, 0, 100) . "\n";

# sample from family_members.annot:
# >up|Q00257|1A12_CUCMA
# CL61	40	424	2.1e-143	pf21ls	1
# >up|Q00379|1A12_CUCPE
# CL61	50	431	1.9e-150	pf21ls	1
# >up|O70212|5HT3A_CAVPO
# PF02931	39	248	3e-73	pf21ls	1
# PF02932	255	481	3e-42	pf21ls	1
# >up|P46098|5HT3A_HUMAN
# PF02931	34	242	5.6e-74	pf21ls	1
# PF02932	249	469	2.5e-46	pf21ls	1

my @entries = split( '>', $domainLocsFileContents);
print STDERR "DEBUG: Found " . ($#entries + 1 - 1) ." entries in $domainLocsFilename\n";

my @seqsMulti;  # taxa names / gis
my $domains = [[[]]];
my @taxaWithRepeatedDomains;
my $numEntriesWithRepeatedDomains = 0;
my @collapsedDAStrs; # indexed by (unique) (multiple) domain architecture (DA) index; values are the (collapsed) domain strings
my $das = [[]];        # first index is (unique) (multiple) domain architecture (DA) index; second indices are the domain strings (so it's a ragged 2D array)
my $uniqueDAIndex = 0;
my $lengthOfHighestUniqueDAIndex = 0;
my $dasBySize = [[]];  # first index is size (+ 2 [because DAs all of multiple members); second index is the $das index
my %uniqueDAs2index;   # keys are unique DA strings (comma seperated list of domains); values are the $das index
my $da2taxa = [[]];  # first index is unique DA index; second index is taxon index (in seqsMulti) (so it's a ragged 2D array)
my @hasRepeats; # index: unique DA index
my @superSetDA; # index: unique DA index
my @numDomainsCollapsed; # index: unique DA index;  number of domains that were collapsed in the seed (or initial) sequence
my @numSequencesWithDomainsCollapsed; # index: unique DA index;  number of member sequences with collapsed domains (including the seed sequence)
my @candidateForQueries;  # index: unique DA index; For choosing queries; 
my @taxon2da; # mapping of taxon to DA
my $numSingleDomains = 0;
my $numSeqsWithMultipleDomainsAllOfTheSameType = 0;

if( $entries[0] !~ m/^\s*$/){
    die( "ERROR: entries[0]: $entries[0]");
}
for ( my $i = 1; $i <= $#entries; $i++){
    if( $entries[$i] !~ m/\n.+\n.+\n/){
	# single-domain
	$numSingleDomains++;
	next;
    }
    my @entry = split( /\s+/, $entries[$i]);
    my $taxon = shift @entry;
    my $potentialTaxonIndex = $#seqsMulti + 1;
    if( ($#entry + 1) % $NUMBER_OF_DOMAIN_FIELDS != 0){
	print STDERR "ERROR: i: $i, mod: ". ($#entry + 1) % $NUMBER_OF_DOMAIN_FIELDS . ", entry: $entries[$i]";
	print STDERR "\t:", join(':', @entry) . ":\n";
	die();
    }

    my $domainsRecorded = 0;
    #$taxaWithRepeatedDomains[ $potentialTaxonIndex] = $FALSE;
    for(my $domainIndex = 0; ($domainIndex + 1) * $NUMBER_OF_DOMAIN_FIELDS <= $#entry + 1; $domainIndex++){
	if( $entry[ $domainIndex * $NUMBER_OF_DOMAIN_FIELDS + $MODE_FIELD] eq "pfam21fs"){
	    # ignore "fs" matches
	    print STDERR "DEBUG: Skipping pfam21fs matches " . join( "\t", subarray(\@entry, $domainIndex * $NUMBER_OF_DOMAIN_FIELDS, $NUMBER_OF_DOMAIN_FIELDS)) ."\n";
	    next;
	}

	for( my $fieldIndex = 0; $fieldIndex < $NUMBER_OF_DOMAIN_FIELDS  &&  $domainIndex * $NUMBER_OF_DOMAIN_FIELDS + $fieldIndex <= $#entry ; $fieldIndex++){
	    $domains->[$potentialTaxonIndex][$domainsRecorded][$fieldIndex] = $entry[ $domainIndex * $NUMBER_OF_DOMAIN_FIELDS + $fieldIndex];
	}
	$domainsRecorded++;
    }

    
    # build DA string
    my $domainIndex = 0;
    my $daStr = $domains->[$potentialTaxonIndex][$domainIndex][$DOMAIN_NAME_FIELD];
    my $daStrFull = $daStr;
    my $numDomainsCollapsed = 0;
    for( $domainIndex = 1; $domainIndex <= $#{$domains->[$potentialTaxonIndex]}; $domainIndex++){
	# collapse consecutive repeated domains
	#print STDERR "Comparing domain names (domainIndex: $domainIndex): $domains->[$potentialTaxonIndex][$domainIndex][$DOMAIN_NAME_FIELD] and $domains->[$potentialTaxonIndex][$domainIndex - 1][$DOMAIN_NAME_FIELD]\n";
	if( ! ($domains->[$potentialTaxonIndex][$domainIndex][$DOMAIN_NAME_FIELD] eq $domains->[$potentialTaxonIndex][$domainIndex - 1][$DOMAIN_NAME_FIELD])){
	    $daStr .= ',' . $domains->[$potentialTaxonIndex][$domainIndex][$DOMAIN_NAME_FIELD];
	}else{
	    #print STDERR "\tcollapsed (numDomainsCollapsed: $numDomainsCollapsed; domainsRecorded: $domainsRecorded)\n";
	    $numDomainsCollapsed++;
	    $domainsRecorded--;
	}
	$daStrFull .= ',' . $domains->[$potentialTaxonIndex][$domainIndex][$DOMAIN_NAME_FIELD];
    }
    
    if( $numDomainsCollapsed == $#{$domains->[$potentialTaxonIndex]} ){
	print STDERR "DEBUG: Skipping because, while there's ".($#{$domains->[$potentialTaxonIndex]} + 1)." domains, they're all the same!\n";
	
	# clear out row in domains multi-array
	$#{$domains->[$potentialTaxonIndex]} = 0;
	
	$numSeqsWithMultipleDomainsAllOfTheSameType++;
	next;
    }
    
    push( @seqsMulti, $taxon);

    if( ! exists( $uniqueDAs2index{$daStr})){
	# new DA
	# check for repeated domains (for information purposes only)
	my $noRepeatedDomains = $TRUE;
	$hasRepeats[ $uniqueDAIndex] = $FALSE;
	$superSetDA[ $uniqueDAIndex] = $FALSE;	
	TOP: for(my $domainIndex = 0; $domainIndex <= $#{$domains->[$#seqsMulti]}; $domainIndex++){
	    for( my $domain2Index = 0; $domain2Index < $domainIndex; $domain2Index++){
		if( $domains->[$#seqsMulti][$domainIndex ][$DOMAIN_NAME_FIELD] eq
		    $domains->[$#seqsMulti][$domain2Index][$DOMAIN_NAME_FIELD] ){
		    $numEntriesWithRepeatedDomains++;
		    $noRepeatedDomains = $FALSE;
		    $hasRepeats[ $uniqueDAIndex] = $TRUE;
		    last TOP;
		}
	    }
	}

	#$candidateForQueries[ $uniqueDAIndex] = $noRepeatedDomains;
	$candidateForQueries[ $uniqueDAIndex] = $TRUE;

	print STDERR "Found unique DA: $daStr (index: $uniqueDAIndex; size: $domainsRecorded; taxon: " . $seqsMulti[$#seqsMulti] ."; multi-domain seq index: $#seqsMulti; noRepeatedDomains: $noRepeatedDomains; domains collapsed: $numDomainsCollapsed)\n";

	if( $domainsRecorded < 0 ){
	    die();
	}

	$uniqueDAs2index{$daStr} = $uniqueDAIndex;
	$collapsedDAStrs[$uniqueDAIndex] = $daStr;
	for( my $domainIndex = 0; $domainIndex <= $#{$domains->[$#seqsMulti]}; $domainIndex++){
	    $das->[$uniqueDAIndex][$domainIndex] = $domains->[$#seqsMulti][$domainIndex][$DOMAIN_NAME_FIELD];
	}

	# record the DA (based on its size)
	push( @{$dasBySize->[$domainsRecorded - 2]}, $uniqueDAIndex);  # "- 2" because we're only considering multi-domains (and

	push( @{$da2taxa->[$uniqueDAIndex]}, $#seqsMulti);  # associate unique DA with taxon
	push( @taxon2da, $uniqueDAIndex);

	$numDomainsCollapsed[ $uniqueDAIndex] = $numDomainsCollapsed;
	if( $numDomainsCollapsed > 0){
	    $numSequencesWithDomainsCollapsed[ $uniqueDAIndex] = 1;  # count the seed sequence
	    if( $DEBUG ){
		print STDERR "DEBUG: COLLAPSED: $numDomainsCollapsed:\n$daStrFull\n$daStr\n";
	    }
	}else{
	    $numSequencesWithDomainsCollapsed[ $uniqueDAIndex] = 0;
	}
	
	$uniqueDAIndex++;
    }else{
	print STDERR "Previously found unique DA: $daStr (index: $uniqueDAs2index{$daStr}; taxon: $seqsMulti[$#seqsMulti]; multi-domain seq index: $#seqsMulti)";
	push( @{$da2taxa->[$uniqueDAs2index{$daStr}]}, $#seqsMulti);  # associate unique DA with taxon
	push( @taxon2da, $uniqueDAs2index{$daStr});
	if( $numDomainsCollapsed != $numDomainsCollapsed[ $uniqueDAs2index{$daStr} ] ){
	    $numSequencesWithDomainsCollapsed[ $uniqueDAs2index{$daStr} ]++;
	    print STDERR " (num collapsed: $numDomainsCollapsed != $numDomainsCollapsed[ $uniqueDAs2index{$daStr} ])";
	}
	print STDERR "\n";
    }
}    

# bucket number of sequences with domains collapsed
my $numDAsWithOnlyTheSeedSequenceCollapsed = 0;
my %numSequencesWithDomainsCollapsedHisto;
for( my $i = 0; $i < $uniqueDAIndex; $i++){
    $numSequencesWithDomainsCollapsedHisto{ $numSequencesWithDomainsCollapsed[ $i] }++;
    
    if( $numSequencesWithDomainsCollapsed[ $i] == 1 &&
	$numDomainsCollapsed[ $i] != 0){
	$numDAsWithOnlyTheSeedSequenceCollapsed++;
    }
}
print STDERR "Before criteria for query DAs has been applied ...\n";
print STDERR "$numDAsWithOnlyTheSeedSequenceCollapsed DAs have only the seed sequence collapsed\n";
print STDERR "Histogram of the number of sequences in a DA with collapsed DAs:\n";
print STDERR "NumSeqs\tOccurrences\n";

foreach my $numSeqs (sort { $a <=> $b } keys %numSequencesWithDomainsCollapsedHisto){
    print STDERR "$numSeqs\t$numSequencesWithDomainsCollapsedHisto{$numSeqs}\n";
}

$lengthOfHighestUniqueDAIndex = length( $uniqueDAIndex);

if( $VERBOSE){
    print STDERR "Summary: Found $numSingleDomains single domain entries\n";
    print STDERR "Summary: Found $numSeqsWithMultipleDomainsAllOfTheSameType sequences with multiple domains, but they're all of the same type\n";
    print STDERR "Summary: Found " . ($#{$domains} + 1) . " multi-domain entries (with at least two types of domains) ($numEntriesWithRepeatedDomains with repeated domains)\n";
    print STDERR "Summary: Found $uniqueDAIndex unique DAs\n";
    #print STDERR "Found the following Taxa: ". join(',', @seqsMulti)."\n";
}

my $defaultSeed = 999219268; # based on a random seed
my $seed = srand( $defaultSeed);
print STDERR "seed: $seed\n";

# 1) Randomly choose a sequence from every base DA
# 2) Find supersets
my @uniqueDAIndex2RandomQuery;
my @uniqueDAIndex2ShortestQuery;
my @uniqueDAIndex2LongestQuery;
my @allQueriesGreaterThanMax; # index: unique DA index
for( my $daIndex = 0; $daIndex <= $#taxon2da; $daIndex++){
    $allQueriesGreaterThanMax[$daIndex] = 2; # unknown
}

my $numOfSuperSets = 0;  # number of times a superset is found (different than numOfSuperSetDAs)
for( my $sizeIndex = 0; $sizeIndex <= $#{$dasBySize}; $sizeIndex++){
    print STDERR "Working on DAs of size " . ($sizeIndex + 2) . " (" . ($#{$dasBySize->[$sizeIndex]} + 1)." unique DAs)\n";

    for( my $daSizeIndex = 0; $daSizeIndex <= $#{$dasBySize->[$sizeIndex]}; $daSizeIndex++){
	print STDERR "\tdaSizeIndex: $daSizeIndex";
	my $uniqueDAIndex = $dasBySize->[$sizeIndex][$daSizeIndex];
	print STDERR "; uniqueDAIndex: $uniqueDAIndex (".( $#{$da2taxa->[$uniqueDAIndex]} + 1)." Taxa)";

	# disqualify DAs with only 1 member
	if( ( $#{$da2taxa->[$uniqueDAIndex]} + 1) <= 1){
	    $candidateForQueries[$uniqueDAIndex] = $FALSE;
	}
	
	# ignore DAs with only 1 member or have a length that's too small or too long
	if( $candidateForQueries[$uniqueDAIndex] == $TRUE){
	    # Make an array of random indices by shuffling an in order one:
	    my @randIndices;
	    for( my $i = 0; $i <= $#{$da2taxa->[$uniqueDAIndex]}; $i++){
		$randIndices[$i] = $i;
	    }
	    fisherYatesShuffle( \@randIndices);

	    my $randomIndex = 0;
	    my $randomTaxonIndex;
	    for( $randomIndex = 0; $randomIndex <= $#randIndices; $randomIndex++){
		$randomTaxonIndex = $da2taxa->[$uniqueDAIndex][$randIndices[$randomIndex]];
		my $randomTaxon = $seqsMulti[$randomTaxonIndex];
		my $randomTaxonLen = getLen( $randomTaxon);

		print STDERR "; randomTaxonIndex: $randomTaxonIndex ($randomTaxon (len: $randomTaxonLen);)";
		if( $randomTaxonLen < $MIN_QUERY_LENGTH){
		    print STDERR "(too short)";
		}elsif( $randomTaxonLen > $MAX_QUERY_LENGTH){
		    print STDERR "(too long)";
		}else{
		    # found match
		    last;
		}
	    }

	    if( $randomIndex > $#randIndices){
		print STDERR "None of the seqs have a length between $MIN_QUERY_LENGTH..$MAX_QUERY_LENGTH!\n";
		$candidateForQueries[$uniqueDAIndex] = $FALSE;
		$allQueriesGreaterThanMax[ $uniqueDAIndex] = $TRUE;
		next;
	    }
	    $allQueriesGreaterThanMax[ $uniqueDAIndex] = $FALSE;
	    $uniqueDAIndex2RandomQuery[$uniqueDAIndex] = $randomTaxonIndex;

	    # find longest and shortest queries
	    my $shortestIndex = 0;
	    my $shortestLen = getLen( $seqsMulti[$da2taxa->[$uniqueDAIndex][$shortestIndex]]);
	    my $longestIndex = 0;
	    my $longestLen = getLen( $seqsMulti[$da2taxa->[$uniqueDAIndex][$longestIndex]]);
	    
	    for( my $taxonIndex = 1; $taxonIndex <= $#{$da2taxa->[$uniqueDAIndex]}; $taxonIndex++){
		my $taxonIndexLen = getLen( $seqsMulti[$da2taxa->[$uniqueDAIndex][$taxonIndex]]);
		if( $taxonIndexLen < $shortestLen){
		    $shortestIndex = $taxonIndex;
		    $shortestLen = $taxonIndexLen;
		}elsif( $taxonIndexLen > $longestLen){
		    $longestIndex = $taxonIndex;
		    $longestLen = $taxonIndexLen;
		}
	    }
	    my $shortestTaxonIndex = $da2taxa->[$uniqueDAIndex][$shortestIndex];
	    my $shortestTaxon = $seqsMulti[$shortestTaxonIndex];
	    my $longestTaxonIndex = $da2taxa->[$uniqueDAIndex][$longestIndex];
	    my $longestTaxon = $seqsMulti[$longestTaxonIndex];
	    
	    print STDERR ", shortestTaxonIndex: $shortestTaxonIndex ($shortestTaxon (len: $shortestLen))";
	    print STDERR ", longestTaxonIndex: $longestTaxonIndex ($longestTaxon (len: $longestLen))";
	    $uniqueDAIndex2ShortestQuery[$uniqueDAIndex] = $shortestTaxonIndex;
	    $uniqueDAIndex2LongestQuery[$uniqueDAIndex] = $longestTaxonIndex;

	    print STDERR ": $collapsedDAStrs[$uniqueDAIndex]\n";
	    
	    # Look through all DAs that are bigger then this DA to see if one
	    my $foundSuperSet = $FALSE;
	    for( my $size2Index = $sizeIndex + 1; $size2Index <= $#{$dasBySize}; $size2Index++){
		#print STDERR "\t\tLooking for DA supersets (size: ". ($size2Index + 2). ")\n";
		for( my $daSize2Index = 0; $daSize2Index <= $#{$dasBySize->[$size2Index]}; $daSize2Index++){
		    my $uniqueDA2Index = $dasBySize->[$size2Index][$daSize2Index];
		    #print STDERR "\t\t\tAnalyzing daSize2Index: $daSize2Index (uniqueDAIndex: $uniqueDA2Index): ".join(',', @{$das->[$uniqueDA2Index]})."\n";
		    my $domain1Index = 0;
		    for( my $domain2Index = 0; $domain2Index < $size2Index + 2; $domain2Index++){
			#print STDERR "\t\t\t\tdomain2Index: $domain2Index\n";
			if( $das->[$uniqueDAIndex][$domain1Index] eq $das->[$uniqueDA2Index][$domain2Index]){
			    if( $domain1Index == $sizeIndex + 2 - 1){
				#print STDERR "\t\t\tFound superset: ". join(',', @{$das->[$uniqueDA2Index]})." (uniqueDA2Index: $uniqueDA2Index)\n";
				for( my $taxonIndex = 0; $taxonIndex <= $#{$da2taxa->[$uniqueDAIndex]}; $taxonIndex++){
				    $taxon2da[$da2taxa->[$uniqueDAIndex][$taxonIndex]] .= ",$uniqueDA2Index";
				}
				push( @{$da2taxa->[$uniqueDAIndex]}, @{$da2taxa->[$uniqueDA2Index]});  # add Taxa from the superset DA
				# $da2taxa->[$uniqueDA2Index] = ();  # a DA may be a superset to more than one DA (base set)
				# splice( @{$dasBySize->[$size2Index]}, $daSize2Index, 1); # a DA may be a superset to more than one DA (base set)

				$numOfSuperSets++;
				$superSetDA[$uniqueDA2Index] = $TRUE;
				#$candidateForQueries[$uniqueDA2Index] = $FALSE;
				$foundSuperSet = $TRUE;
				last;
			    }
			    $domain1Index++;
			}
		    }
		}
	    }
	    if( $foundSuperSet){
		print STDERR "\t\tNow with ".( $#{$da2taxa->[$uniqueDAIndex]} + 1)." Taxa\n";
	    }

	}else{
	    # Determine if all of the members of the DA are too short or too long (just stats)
	    my $i = 0;
	    for( ; $i <= $#{$da2taxa->[$uniqueDAIndex]}; $i++){
		my $taxonIndex = $da2taxa->[$uniqueDAIndex][$i];
		my $taxon = $seqsMulti[$taxonIndex];
		my $taxonLen = getLen( $taxon);

		print STDERR "; taxonIndex: $taxonIndex ($taxon (len: $taxonLen);)";
		if( $taxonLen < $MIN_QUERY_LENGTH){
		    print STDERR "(too short)";
		}elsif( $taxonLen > $MAX_QUERY_LENGTH){
		    print STDERR "(too long)";
		}else{
		    # found match
		    last;
		}
	    }

	    if( $i > $#{$da2taxa->[$uniqueDAIndex]}){
		print STDERR "None of the seqs have a length between $MIN_QUERY_LENGTH..$MAX_QUERY_LENGTH!\n";
		$allQueriesGreaterThanMax[ $uniqueDAIndex] = $TRUE;
	    }else{
		$allQueriesGreaterThanMax[ $uniqueDAIndex] = $FALSE;
	    }
	    
	    print STDERR "\n";
	}
    }
}
#
# Stats
#

# indices:
# 1) notMoreThan1Member (0: good (potential DA for query), 1: singleton DA)
# 2) allQueriesGreaterThanMax (0: good (potential DA for query), 1: bad)
# 3) hasRepeats (0, 1)
# 4) isSuperSet (0, 1)
my $statsSummary = [[[[]]]]; 
for( my $notMoreThan1Member = 0; $notMoreThan1Member <= 1; $notMoreThan1Member++){
    for( my $allQueriesGreaterThanMax = 0; $allQueriesGreaterThanMax <= 1; $allQueriesGreaterThanMax++){
	for( my $isSuperSet = 0; $isSuperSet <= 1; $isSuperSet++){
	    for( my $isRepeat = 0; $isRepeat <= 1; $isRepeat++){
		$statsSummary->[$notMoreThan1Member][$allQueriesGreaterThanMax][$isSuperSet][$isRepeat] = 0;
	    }
	}
    }
}

my $numOfDAsWithSingleTaxon = 0;
my $numOfDAsWithSingleTaxonSuperSet = 0;
my $numOfDAsWithSingleTaxonWithRepeatDomains = 0;
my $numOfDAsWithSingleTaxonRepeatSuperSet = 0;
my $numOfSuperSetDAs = 0; # number of DAs that are an ordered superset of 1 or more DAs
my $numOfDAsWithRepeatedDomains = 0;
my $numRepeatSuperSetDAs = 0;
for( my $daIndex = 0; $daIndex <= $#candidateForQueries; $daIndex++){
    my $notMoreThan1Member = 0;  # 0: Has more than 1 member; 1: Singleton DA
    if( ( $#{$da2taxa->[$daIndex]} + 1) <= 1){
	$numOfDAsWithSingleTaxon++;
	if( $superSetDA[$daIndex]){
	    $numOfDAsWithSingleTaxonSuperSet++;
	}
	if($hasRepeats[$daIndex]){
	    $numOfDAsWithSingleTaxonWithRepeatDomains++;
	    if( $superSetDA[$daIndex] ){
		$numOfDAsWithSingleTaxonRepeatSuperSet++;
	$numRepeatSuperSetDAs++;
	    }
	}
	$notMoreThan1Member = 1;
    }

    $numOfSuperSetDAs += $superSetDA[$daIndex];
    
    $numOfDAsWithRepeatedDomains += $hasRepeats[$daIndex];

    if( $allQueriesGreaterThanMax[$daIndex] == 2 ){
	die("ERROR: DA index $daIndex is unknown for allQueriesGreaterThanMax!");
    }
    
    $statsSummary->[$notMoreThan1Member][$allQueriesGreaterThanMax[$daIndex]][$superSetDA[$daIndex]][$hasRepeats[$daIndex]]++;

    # if( $candidateForQueries[ $daIndex] != (($notMoreThan1Member == 0) && ($allQueriesGreaterThanMax[$daIndex] == 0)) && ($superSetDA[$daIndex] == $FALSE) ){
    # 	die( "ERROR: candidateForQueries[ $daIndex] ($candidateForQueries[ $daIndex]) != (($notMoreThan1Member == 0) && ($allQueriesGreaterThanMax[$daIndex] == 0)) && ($superSetDA[$daIndex] == $FALSE) )");
    # }
    if( $candidateForQueries[ $daIndex] != (($notMoreThan1Member == 0) && ($allQueriesGreaterThanMax[$daIndex] == 0)) ){
	die( "ERROR: candidateForQueries[ $daIndex] ($candidateForQueries[ $daIndex]) != (($notMoreThan1Member == 0) && ($allQueriesGreaterThanMax[$daIndex] == 0)) )");
    }
}
print STDERR "DA STAT SUMMARY: #notMoreThan1Member,allQueriesGreaterThanMax,isSuperSet,hasRepeats,count\n";
for( my $notMoreThan1Member = 0; $notMoreThan1Member <= 1; $notMoreThan1Member++){
    for( my $allQueriesGreaterThanMax = 0; $allQueriesGreaterThanMax <= 1; $allQueriesGreaterThanMax++){
	for( my $isSuperSet = 0; $isSuperSet <= 1; $isSuperSet++){
	    for( my $isRepeat = 0; $isRepeat <= 1; $isRepeat++){
		print STDERR "DA STAT SUMMARY: $notMoreThan1Member,$allQueriesGreaterThanMax,$isSuperSet,$isRepeat,$statsSummary->[$notMoreThan1Member][$allQueriesGreaterThanMax][$isSuperSet][$isRepeat]\n";
	    }
	}
    }
}

print STDERR "Number of unique DAs: ".($#candidateForQueries + 1)."\n";
print STDERR "Number of DAs with only one sequence: $numOfDAsWithSingleTaxon ($numOfDAsWithSingleTaxonSuperSet also are a super set; $numOfDAsWithSingleTaxonWithRepeatDomains also have a repeated domain; $numOfDAsWithSingleTaxonRepeatSuperSet are both a super set and have a repeated domain)\n";
print STDERR "Number of DAs with repeats: $numOfDAsWithRepeatedDomains (out of ".($#hasRepeats + 1).")\n";
print STDERR "Number of DAs that are ordered supersets: $numOfSuperSetDAs (out of ".($#superSetDA + 1).") (number of times a DA is a superset: $numOfSuperSets [a DA can have multiple ordered superset DAs])\n";
print STDERR "Number of ordered superset DAs with repeats: $numRepeatSuperSetDAs\n";

# bucket number of sequences with domains collapsed
$numDAsWithOnlyTheSeedSequenceCollapsed = 0;
%numSequencesWithDomainsCollapsedHisto = ();
for( my $daIndex = 0; $daIndex <= $#candidateForQueries; $daIndex++){
    if( $candidateForQueries[ $daIndex] ){  # only count DAs that are candidates for queries
	$numSequencesWithDomainsCollapsedHisto{ $numSequencesWithDomainsCollapsed[ $daIndex] }++;
	
	if( $numSequencesWithDomainsCollapsed[ $daIndex] == 1 &&
	    $numDomainsCollapsed[ $daIndex] != 0){
	    $numDAsWithOnlyTheSeedSequenceCollapsed++;
	}
    }
}
print STDERR "AFTER criteria for query DAs has been applied ...\n";
print STDERR "$numDAsWithOnlyTheSeedSequenceCollapsed DAs have only the seed sequence collapsed\n";
print STDERR "Histogram of the number of sequences in a DA with collapsed DAs:\n";
print STDERR "NumSeqs\tOccurrences\n";

foreach my $numSeqs (sort { $a <=> $b } keys %numSequencesWithDomainsCollapsedHisto){
    print STDERR "$numSeqs\t$numSequencesWithDomainsCollapsedHisto{$numSeqs}\n";
}


print STDERR "Purge queries that don't have other sequences with the same DA and print relevance info file\n";

# save information about each multiple DA to a file
open( DA_INFO, ">$daInformationFilename") or die( "ERROR: Can not open / create daInformationFilename, \"$daInformationFilename\": $!");
print DA_INFO "#DA\tDomains\tCollapsed_Domains\tTaxa\n";
for( my $daIndex = 0; $daIndex <= $#{$da2taxa}; $daIndex++){
    my $daId = "da" . padNumber( $daIndex, $lengthOfHighestUniqueDAIndex + 1);
    my $taxonIndex = 0;
    my $taxonStr = $seqsMulti[ $da2taxa->[$daIndex][$taxonIndex]];
    $taxonIndex++;
    for( ; $taxonIndex <= $#{$da2taxa->[$daIndex]}; $taxonIndex++){
	$taxonStr .= ",$seqsMulti[ $da2taxa->[$daIndex][$taxonIndex]]";
    }
    my $entry = "$daId\t$collapsedDAStrs[$daIndex]\t";
    if( $numSequencesWithDomainsCollapsed[ $daIndex] > 0){
	$entry .= "T";
    }else{
	$entry .= "F";
    }
    $entry .= "\t$taxonStr\n";
    print DA_INFO $entry;
}
close( DA_INFO);


open( RELEVANCE, ">$relevanceFilename") or die( "ERROR: Can not open / create relevanceFilename, \"$relevanceFilename\": $!");
print RELEVANCE "#Taxon\tDA\tSuperset_DAs\n";

for( my $taxonIndex = 0; $taxonIndex <= $#taxon2da; $taxonIndex++){
    my @allDas = split( /,/, $taxon2da[$taxonIndex]);
    my $daId = "da" . padNumber( $allDas[0], $lengthOfHighestUniqueDAIndex + 1);
    my $relStr = "$seqsMulti[$taxonIndex]\t$daId";
    
    my %uniqueList;
    for( my $i = 1; $i <= $#allDas; $i++){
	$uniqueList{ $allDas[$i]} = 1;
    }
    my $seperator = "\t";
    foreach my $daIndex (sort { $a <=> $b} keys %uniqueList){
	my $daId = "da" . padNumber( $daIndex, $lengthOfHighestUniqueDAIndex + 1);
	$relStr .= "$seperator$daId";
	$seperator = ',';
    }
    print RELEVANCE "$relStr\n";
}
close( RELEVANCE);
    
open( QUERIES, ">$queriesFilename") or die( "ERROR: Can not open / create queriesForDAs \"$queriesFilename\": $!");
print QUERIES "#DA_ID\tRandom_Taxon\tShortest_Taxon\tLongest_Taxon\n";
# only need to go until the last DA that is a candidate for a query (and not all the way to $#candidateForQueries)
for( my $daIndex = 0; $daIndex <= $#uniqueDAIndex2RandomQuery; $daIndex++){
    #print STDERR "uniqueDAIndex2RandomQuery: daIndex $daIndex: ";
    
    if( $candidateForQueries[$daIndex]){
	my $daId = "da" . padNumber( $daIndex, $lengthOfHighestUniqueDAIndex + 1);	
	print QUERIES "$daId\t$seqsMulti[$uniqueDAIndex2RandomQuery[$daIndex]]\t$seqsMulti[$uniqueDAIndex2ShortestQuery[$daIndex]]\t$seqsMulti[$uniqueDAIndex2LongestQuery[$daIndex]]\n";
    }else{
	print STDERR "Ignoring daIndex $daIndex: ";
	if( defined( $hasRepeats[$daIndex])){
	    print STDERR "hasRepeats: $hasRepeats[$daIndex]; ";
	}
	if( defined($superSetDA[$daIndex])){
	    print STDERR "superSetDA: $superSetDA[$daIndex];";
	}
	print STDERR "(".($#{$da2taxa->[$daIndex]} + 1)." Taxa)\n";
    }
}
close( QUERIES);

print STDERR "Done ($0)\n";
exit(0);


sub subarray{
    my ($arrayRef, $startIndex, $length) = @_;

    my @retval;
    for( my $i = $startIndex; $i < $startIndex + $length && $i <= $#{$arrayRef}; $i++){
	push( @retval, $arrayRef->[$i]);
    }
    return @retval;
}

sub padNumber{
    my ($num, $digits) = @_;
    my $retval = "";
    for( my $len = length( $num); $len < $digits; $len++){
	$retval .= "0";
    }
    return $retval . $num;
}


sub getLen{
    my ($taxon) = @_;
    if( ! exists( $seqLens{$taxon})){
	die("\nERROR: Could not find hash entry for $taxon!");
    }
    return $seqLens{$taxon};
}


__END__

=head1 NAME

getDAs.pl - Given a FASTA file and the location of domains, determine all of the
            Domain Architectures (DAs) present and assigns each sequence to a
            DA.  Additionally, potential queries are filtered (e.g., by length).

=head1 SYNOPSIS

getDAs.pl [INPUTS] [OUTPUTS] [OPTIONS]

 Inputs (all are optional):
   --fasta <FASTA formated filename (defaults to library_all_domains.fa)>
   --domainLocs <filename with locations of domains (defaults to
                 family_members-nonOverlapping.annot)>

 Outputs (all are optional):
   --da <filename of TAB-delimited DA IDs and associated taxa (along with
         domains and the if there are collapsed domains) (defaults to
         daInfo.tab)>
   --rel <filename of TAB-delimited Taxon and DA (relevance) information
         (defaults to taxon2da.tab)>
   --queries <filename of taxon (random, shortest and longest) per each DA in a
              TAB-delimited formated file (defaults to queriesPerDA.tab)>
    
 Options (all are optional):
   -v | --verbose  Display verbose output
   -d | --debug  Display debugging information
   --help  Brief help message
   --man  Full documentation

=cut
