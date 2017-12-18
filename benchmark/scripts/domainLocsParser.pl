#!/usr/bin/env perl

# Description: Creates domainQueries_multi.tab and domainQueries_single.tab  (or 2nd & 3rd command-line args) from family_members-nums.annot (or 1st command-line arg)

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

BEGIN
{
    # Look for perl modules in ./perl
    unshift( @INC, "./perl"); 
}

use HDC::Common qw(:ALL);


#my $NUMBER_OF_FIELDS = 7;
my $NUMBER_OF_FIELDS = 6;
my $NAME_FIELD = 0;        # Assumes the domain name is the first field
my $START_FIELD = 1;         # Assumes the domain starting position is the second field
my $LAST_FIELD = 2;          # Assumes the domain ending position is the third field
my $UNKNOWN_FIELD = 3;
my $MODE_FIELD = 4;               # Assumes the mode is the fifth field
#my $HAS_REPEATS_COL = $MODE_FIELD + 1;
#my $OVERLAP_TYPE_COL = $HAS_REPEATS_COL + 1;

#
# Command line options
#

# Configuration defaults, can be overriden using command line arguments
# $DEBUG = 1;
my $domainLocsFilenameIn = "family_members.annot";
my $domainLocsFilenameOut = "";
my $singleDomainOutputFilename = "domainQueries_single.tab";
my $multiDomainOutputFilename = "domainQueries_multi.tab";
my $removedTaxaFilename;
my $man = 0;
my $help = 0;

my $usageStr = "Usage: $0\n".
    "-in=<original domain locations filename> (defaults to $domainLocsFilenameIn)\n".
    "[-out=<domain locations output filename> (defaults to STDOUT)]\n".
    "[-single=<filename> (Single-domain info output) (defaults to $singleDomainOutputFilename)\\n".
    "[-multi=<filename> (Multi-domain info output) (defaults to $multiDomainOutputFilename)\\n".
    "[-removed=<filename> (List of taxa removed (b/c they have domains that overlap))]\n".
    "-verbose\n".
    "-debug\n".
    "-help";


unless (GetOptions("verbose+"            => \$VERBOSE,
		   "debug+"              => \$DEBUG,
		   "in=s"                => \$domainLocsFilenameIn,
		   "out=s"               => \$domainLocsFilenameOut,
		   "single=s"            => \$singleDomainOutputFilename,
		   "multi=s"             => \$multiDomainOutputFilename,
		   "removed=s"           => \$removedTaxaFilename,
		   'help|?'              => \$help,
		   'man'                 => \$man)){
    print STDERR "$usageStr\n";
    exit (0);
}

VERBOSE("domainLocsFilenameIn: $domainLocsFilenameIn");
VERBOSE("domainLocsFilenameOut: $domainLocsFilenameOut");
VERBOSE("singleDomainOutputFilename: $singleDomainOutputFilename");
VERBOSE("multiDomainOutputFilename: $multiDomainOutputFilename");
VERBOSE("removedTaxaFilename: $removedTaxaFilename");
VERBOSE("VERBOSE: $VERBOSE");
VERBOSE("DEBUG: $DEBUG");

my $domainLocsFileContents = fileToString( $domainLocsFilenameIn);
#print "domainLocsFileContents:\n" . substr( $domainLocsFileContents, 0, 100) . "\n";

my @entries = split( '>', $domainLocsFileContents);
my @seqsMulti;   # taxa names / gis
my @seqsSingle;  # taxa names / gis
my $domains = [[[]]];
my $singleDomainLocsStr = "";
my $multiDomainLocsStr = "";
for ( my $i = 0; $i <= $#entries; $i++){
    my @entry = split( /\s+/, $entries[$i]);
    my $gi = shift @entry;
    if( $entries[$i] =~ m/\n.+\n.+\n/){
	# multi-domain
	push( @seqsMulti, $gi);
	if( ($#entry + 1) % $NUMBER_OF_FIELDS != 0){
	    print STDERR "ERROR: i: $i, mod: ". ($#entry + 1) % $NUMBER_OF_FIELDS . ", entry: $entries[$i]";
	    print STDERR "\t:", join(':', @entry) . ":\n";
	    die();
	}

	my $domainsRecorded = 0;
	for(my $domainIndex = 0; ($domainIndex + 1) * $NUMBER_OF_FIELDS <= $#entry + 1; $domainIndex++){
	    # #if( $entry[ $domainIndex * $NUMBER_OF_FIELDS + $MODE_FIELD] eq "pfam21fs"){
	    # if( $entry[ $domainIndex * $NUMBER_OF_FIELDS + $MODE_FIELD] eq "pf21fs"){
	    # 	# ignore "fs" matches
	    # 	print STDERR "DEBUG: Skipping " . join( "\t", subarray(\@entry, $domainIndex * $NUMBER_OF_FIELDS, $NUMBER_OF_FIELDS)) ."\n";
	    # 	next;
	    # }
	    for( my $fieldIndex = 0; $fieldIndex < $NUMBER_OF_FIELDS  &&  $domainIndex * $NUMBER_OF_FIELDS + $fieldIndex <= $#entry ; $fieldIndex++){
		$domains->[$#seqsMulti][$domainsRecorded][$fieldIndex] = $entry[ $domainIndex * $NUMBER_OF_FIELDS + $fieldIndex];
	    }
	    
	    $domainsRecorded++;
	}

    }else{
	if( $#entry > 0){
	    push( @seqsSingle, $gi);
	    $singleDomainLocsStr .= ">$entries[$i]";
	}
    }
}    

if( $VERBOSE){
    print STDERR "Found " . ($#{$domains} + 1) . " multi-domain entries (before addressing overlapping domains)\n";
}

my @taxaRemoved;

open( MULTI, ">$multiDomainOutputFilename") or die( "ERROR: Can not open / create multiDomainOutputFilename, \"$multiDomainOutputFilename\": $!");
print MULTI "#GI\tNumber_of_domains\tRepeated_domain?\tEmbedded_domain(s)?\n";
SEQ: for( my $seqIndex = 0; $seqIndex <= $#{$domains}; $seqIndex++){
    #DEBUG("Processing multi-domain seq $seqsMulti[$seqIndex] (seqIndex: $seqIndex). . .");
    my $hasRepeatedDomain = $FALSE;
    my $overlapType = "";
    my $numEmbeddedDomains = 0;
    for( my $domain1Index = 0; $domain1Index <= $#{$domains->[$seqIndex]}; $domain1Index++){
	for( my $domain2Index = 0; $domain2Index < $domain1Index && $domain2Index < $#{$domains->[$seqIndex]}; $domain2Index++){  # 2nd continuation criterion necessary if we spliced out the last one
	    if( $domains->[$seqIndex][$domain1Index][$NAME_FIELD] eq $domains->[$seqIndex][$domain2Index][$NAME_FIELD] ){
		$hasRepeatedDomain = $TRUE;
		#last;
	    }
	    # test if this domain overlaps with any of the previous domains
	    #
	    # 5 Potential overlap scenarios between domain 1 & domain 2:
	    #
	    # Domain1:      +---+
	    # Domain2A: +-+            No overlap
	    # Domain2B:    +-+         Overlap (side)
	    # Domain2C:      +-+       Overlap (sub)
	    # Domain2D:        +-+     Overlap (side)
	    # Domain2E:           +-+  No overlap
	    # Domain2F:   +-------+    Overlap (sub)
	    if( $domains->[$seqIndex][$domain1Index][$START_FIELD] <= $domains->[$seqIndex][$domain2Index][$LAST_FIELD] &&
		$domains->[$seqIndex][$domain1Index][$LAST_FIELD]  >= $domains->[$seqIndex][$domain2Index][$START_FIELD] ){
		# has overlapping domain

		my $startPos = MAX( $domains->[$seqIndex][$domain1Index][$START_FIELD], $domains->[$seqIndex][$domain2Index][$START_FIELD]);
		my $lastPos  = MIN( $domains->[$seqIndex][$domain1Index][$LAST_FIELD],  $domains->[$seqIndex][$domain2Index][$LAST_FIELD]);
		
		my $domainOverlap = $lastPos - $startPos + 1;
		my $domain1Len = $domains->[$seqIndex][$domain1Index][$LAST_FIELD] - $domains->[$seqIndex][$domain1Index][$START_FIELD] + 1;
		my $domain2Len = $domains->[$seqIndex][$domain2Index][$LAST_FIELD] - $domains->[$seqIndex][$domain2Index][$START_FIELD] + 1;
		my $overlapMin = sprintf("%.0f", $domainOverlap / MAX( $domain1Len, $domain2Len) * 100.0);
		my $overlapMax = sprintf("%.0f", $domainOverlap / MIN( $domain1Len, $domain2Len) * 100.0);
		
		# determine type of overlap
		if( ($domains->[$seqIndex][$domain1Index][$START_FIELD] > $domains->[$seqIndex][$domain2Index][$START_FIELD] &&
		     $domains->[$seqIndex][$domain1Index][$LAST_FIELD]  > $domains->[$seqIndex][$domain2Index][$LAST_FIELD] ) ||
		    ($domains->[$seqIndex][$domain1Index][$START_FIELD] < $domains->[$seqIndex][$domain2Index][$START_FIELD] &&
		     $domains->[$seqIndex][$domain1Index][$LAST_FIELD]  < $domains->[$seqIndex][$domain2Index][$LAST_FIELD] )){
		    $overlapType = "side";
		}else{
		    $overlapType = "sub";
		    if( $domains->[$seqIndex][$domain1Index][$NAME_FIELD] eq $domains->[$seqIndex][$domain2Index][$NAME_FIELD] ){
			$overlapType = "dup";
		    }
		    if( ($domains->[$seqIndex][$domain1Index][$MODE_FIELD] eq "pf21fs" &&
			 $domains->[$seqIndex][$domain2Index][$MODE_FIELD] eq "pf21ls") ||
			($domains->[$seqIndex][$domain1Index][$MODE_FIELD] eq "pf21ls" &&
			 $domains->[$seqIndex][$domain2Index][$MODE_FIELD] eq "pf21fs")){
			$overlapType .= " (fs in ls)";
		    }
		}

		DEBUG( "Overlapping domains in $seqsMulti[$seqIndex] (type: $overlapType):\t$domains->[$seqIndex][$domain1Index][$NAME_FIELD] ($domains->[$seqIndex][$domain1Index][$MODE_FIELD]) $domains->[$seqIndex][$domain1Index][$START_FIELD]..$domains->[$seqIndex][$domain1Index][$LAST_FIELD] ($domain1Len) and $domains->[$seqIndex][$domain2Index][$NAME_FIELD] ($domains->[$seqIndex][$domain2Index][$MODE_FIELD]) $domains->[$seqIndex][$domain2Index][$START_FIELD]..$domains->[$seqIndex][$domain2Index][$LAST_FIELD] ($domain2Len); Overlap $domainOverlap ($overlapMin..$overlapMax\%)");

		if( $overlapType =~ m/^dup/){
		    $numEmbeddedDomains++;
		    DEBUG( "Removing the embedded (smaller) domain from $seqsMulti[$seqIndex] ...");
		    if( $domain1Len < $domain2Len){
			# domain 1 is smaller (embedded)
			splice( $domains->[$seqIndex], $domain1Index, 1);
			$domain1Index--;
			last;
		    }else{
			# domain 2 is smaller (embedded)
			splice( $domains->[$seqIndex], $domain2Index, 1);
			$domain2Index--;
			next;
		    }
		}else{
		    DEBUG("Removing taxon $seqsMulti[$seqIndex] because it has overlapping domains (of type $overlapType)\n");
		    push( @taxaRemoved,  $seqsMulti[$seqIndex]);
		    splice( $domains, $seqIndex, 1);
		    splice( @seqsMulti, $seqIndex, 1);
		    $seqIndex--;
		    next SEQ;
		}
	    }
	}
    }
    if( $overlapType ne ""){
	DEBUG( "");
    }
    if( $#{$domains->[$seqIndex]} + 1 <= 1 ){
	DEBUG("It appears that one or more domains have been removed from $seqsMulti[$seqIndex], leaving it with just ".($#{$domains->[$seqIndex]} + 1) ." domain(s)");
	if( $#{$domains->[$seqIndex]} + 1 == 1 ){
	    push( @seqsSingle, $seqsMulti[$seqIndex]);
	    $singleDomainLocsStr .= ">$seqsMulti[$seqIndex]\n";
	    $singleDomainLocsStr .= join( "\t", @{$domains->[$seqIndex][0]}) . "\n";  # 0 b/c just 1 domain
	}
	splice( $domains, $seqIndex, 1);
	splice( @seqsMulti, $seqIndex, 1);
	$seqIndex--;
	next;	
    }
    print MULTI "$seqsMulti[$seqIndex]\t".($#{$domains->[$seqIndex]} + 1)."\t$hasRepeatedDomain\t$numEmbeddedDomains\n";
    $multiDomainLocsStr .= ">$seqsMulti[$seqIndex]\n";
    for( my $domain1Index = 0; $domain1Index <= $#{$domains->[$seqIndex]}; $domain1Index++){
	$multiDomainLocsStr .= join( "\t", @{$domains->[$seqIndex][$domain1Index]}) . "\n";
    }
}
close( MULTI);

if( $VERBOSE){
    print STDERR "Found " . ($#{$domains} + 1) . " multi-domain entries (after addressing overlapping domains) (and ".($#taxaRemoved + 1). " to remove)\n";
}


if( $VERBOSE){
    print STDERR "Found " . ($#seqsSingle + 1). " single-domain entries\n";
}

open( SINGLE, ">$singleDomainOutputFilename") or die( "ERROR: Can not open / create singleDomainOutputFilename, \"$singleDomainOutputFilename\": $!");
print SINGLE "#GI\n";
for( my $seqIndex = 0; $seqIndex <= $#seqsSingle; $seqIndex++){
    print SINGLE "$seqsSingle[$seqIndex]\n";
}
close( SINGLE);


if( $domainLocsFilenameOut ne ""){
    open( OUT, ">$domainLocsFilenameOut") or die("ERROR: Can not open / create domainLocsFilenameOut \"$domainLocsFilenameOut\": $!");
}else{
    *OUT = *STDOUT;
}

print OUT "$multiDomainLocsStr$singleDomainLocsStr";

if( $domainLocsFilenameOut ne ""){
    close( OUT);
}


if( defined( $removedTaxaFilename) && $removedTaxaFilename ne ""){
    open( TAXA_REMOVED, ">$removedTaxaFilename") or die("ERROR: Can not open / create removed taxa file \"$removedTaxaFilename\": $!");
    print TAXA_REMOVED join("\n", @taxaRemoved) . "\n";
    close(TAXA_REMOVED);
}


VERBOSE("Done ($0)");


sub subarray{
    my ($arrayRef, $startIndex, $length) = @_;

    my @retval;
    for( my $i = $startIndex; $i < $startIndex + $length && $i <= $#{$arrayRef}; $i++){
	push( @retval, $arrayRef->[$i]);
    }
    return @retval;
}

exit(0);
