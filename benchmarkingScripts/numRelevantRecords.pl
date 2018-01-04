#!/usr/bin/perl -w

################################################################################
# File: numRelevantRecords.pl
# Input: Taxon label, relevancy information file.
# Output: Outputs the total number of relevant records "TPs".
################################################################################

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

# NOTE: The term "Index" in a variable name identifies the location (index value: 0 .... len - 1) of the start, end, etc. of a string or array
# NOTE: the terms "superfamily" and "fold" can be substituted with "more strigent criterion" and "less strigent criterion"

#
# To-Do
#
#

use strict;
# use Bio::SearchIO;
# use Bio::SeqIO;

my $FALSE = 0;
my $TRUE = 1;

# types of hits output files:
my $UNKNOWN_HITS_OUTPUT_TYPE = -1;
my $BLAST_HITS_OUTPUT_TYPE = $UNKNOWN_HITS_OUTPUT_TYPE + 1;
my $SIMPLE_HITS_OUTPUT_TYPE = $BLAST_HITS_OUTPUT_TYPE + 1;

my $NEVER_DO_FULL_FILES = $TRUE;  # FALSE means account for multiple HSP per database sequence


my $usage = "Usage: $0:
  [-taxon=<taxon name>] (e.g., d12asa_)
  [-family=<superfamily name>] (e.g., 46458 or cd00001)
  [-fold=<fold name>] (e.g., 46457)
  -rel=<relevance info file (e.g., id2superfamily.tab or id2superfamily_fold.tab)>
  [-usefold] (use fold (secondary) classification information in -rel for relevance; Default: unset) 
  [-nofold] (Ignore fold (secondary) classification information in -rel; Default: set)
  [-ignoreSelfHit] (Ignore database hits matching -taxon; Default: set if -taxon set)
  [-allowSelfHit] (Allow for all possible database hits; Default: unset unless -taxon set)
  [-v] (verbose; Default: unset)
  [-d] (debug; Default: unset)";
#  [-out=<output file>]
#  -total=<total number of sequences in the DB>


my %args;  # keys are: 

for( my $argvIndex = 0; $argvIndex <= $#ARGV; $argvIndex++){
    if( $ARGV[ $argvIndex] =~ m/^-([^\s=]+)=(.*)/){
	$args{lc($1)} = $2; # lc = lowercase
	#print "$0 found arg $1 (=$args{$1})\n";
	splice @ARGV, $argvIndex, 1;
	$argvIndex--;
    }elsif( $ARGV[ $argvIndex] =~ m/^-([^\s]+)/){
	$args{lc($1)} = $FALSE; # lc = lowercase
	#print "$0 found arg $1 (=$args{$1})\n";
	splice @ARGV, $argvIndex, 1;
	$argvIndex--;
    }
}

my $VERBOSE = $FALSE;
my $DEBUG = $FALSE;
my $taxon;
my $superfamily;
my $fold;
my $relevanceInfoFileName;
#my $outputFileName = "";
#my $totalSeqs = -1;
my $useFold = $FALSE;
my $ignoreSelfHit = $FALSE;   # defaults to $TRUE later if $taxon is defined

if( exists($args{"v"}) ||
    exists($args{"verbose"})){
    $VERBOSE = $TRUE;
}

if( exists($args{"d"}) ||
    exists($args{"debug"})){
    $DEBUG = $TRUE;
    $VERBOSE = $TRUE;
}

if( defined($args{"taxon"})){
    $taxon = uc($args{"taxon"});
    if( $DEBUG){ print STDERR "DEBUG: taxon: $taxon\n"; }
}
if( defined($args{"family"})){
    $superfamily = uc($args{"family"});
    if( $DEBUG){ print STDERR "DEBUG: superfamily: $superfamily\n"; }
}
if( defined($args{"fold"})){
    $fold = uc($args{"fold"});
    if( $DEBUG){ print STDERR "DEBUG: fold: $fold\n"; }
}
if( ! defined( $taxon) &&
    ! defined( $superfamily) &&
    ! defined( $fold)){
    print STDERR "\nERROR: Either -taxon, -family or -fold is required!\n\n";
    die( "$usage\n");
}


if( defined( $args{"rel"})){
    $relevanceInfoFileName = $args{"rel"};
}
if( ! defined( $relevanceInfoFileName) ||
    $relevanceInfoFileName eq "" ||
    ! -s $relevanceInfoFileName ){
    print STDERR "ERROR: Bad family,superfamily[,fold] info file \"$relevanceInfoFileName\"!\n";
    die( "$usage\n");
}

#if( defined( $args{"out"})){
#    $outputFileName = $args{"out"};
#}


if( exists( $args{"usefold"})){
    $useFold = $TRUE;
}
if( exists( $args{"nofold"})){
    $useFold = $FALSE;
}

my $criterionStr;  # the criterion used to base relevancy (e.g., the specific superfamily that the taxon is a member of)
if( $useFold){
    $criterionStr = $fold;
}else{
    $criterionStr = $superfamily;
}

if( ! defined( $criterionStr) &&
    ! defined( $taxon)){
    print STDERR "\n\nERROR: criterion for secondary classification information not specificed (useFold: $useFold)!\n\n"; 
    die("$usage\n");
}

if( defined( $taxon)){
    $ignoreSelfHit = $TRUE;
}

if( defined( $args{"allowselfhit"})){
    $ignoreSelfHit = $FALSE;
}elsif( defined( $args{"ignoreselfhit"})){
    $ignoreSelfHit = $TRUE;
    if(	! defined( $taxon)){
	print STDERR "\n\nERROR: -taxon must be specified to ignore self hits\n";
	die( "$usage\n");
    }
}


#
# Read in relevant family members
#

my ($relevantTaxaRef, $tpsCounter, $totalSeqs) = parseTaxaInfo($relevanceInfoFileName, $taxon, $criterionStr, $ignoreSelfHit);
if( $DEBUG){ print STDERR "HyDEBUG: tpsCounter: $tpsCounter (ignoreSelfHit: $ignoreSelfHit)\n"; }

print "$tpsCounter\n";

exit(0);


sub parseTaxaInfo{
    my ($relevanceInfoFileName, $taxon, $criterionStr, $ignoreSelfHit) = @_;

    open( TAXON2SUPERFAMILY, "$relevanceInfoFileName") or die("ERROR: relevanceInfoFileName \"$relevanceInfoFileName\";");
    if( defined($_ = <TAXON2SUPERFAMILY>) &&
	$_ !~ m/^#/){
	die( "ERROR: relevanceInfoFileName \"$relevanceInfoFileName\" missing header;");
    }

    my $taxon2criterionRef = {};  # initialize

    if( $useFold){
	while( <TAXON2SUPERFAMILY>){
	    if( $_ =~ m/^(\S+)\s+(\S+)\s+(\S+)/){
		# has taxon, superfamily and fold information
		#if( exists($taxon2criterionRef->{uc($1)})){
		#    print STDERR "ALERT: Ignoring previous fold ($taxon2criterionRef->{uc($1)}) for $1 in favor of $3!\n";
		#}
		$taxon2criterionRef->{uc($1)} .= uc($3) . ",";
	    }else{
		chomp;
		die( "ERROR: Mal-formed taxon2fold line: \"$_\" in $relevanceInfoFileName;");
	    }
	}
    }else{
	while( <TAXON2SUPERFAMILY>){
	    if( $_ =~ m/^(\S+)\s+(\S+)/){
		# has taxon and superfamily information
		#if( exists($taxon2criterionRef->{uc($1)})){
		#    print STDERR "ALERT: Ignoring previous superfamily ($taxon2criterionRef->{uc($1)}) for $1 in favor of $2!\n";
		#}
		$taxon2criterionRef->{uc($1)} .= uc($2) . ",";
	    }else{
		chomp;
		die( "ERROR: Mal-formed taxon2superfamily line: \"$_\" in $relevanceInfoFileName;");
	    }
	}
    }

    if( ! defined( $criterionStr)){
	# criterion wasn't set on the command-line
	# look up the criterion from the parsed information
	if( ! defined( $taxon2criterionRef->{ $taxon})){
	    die( "ERROR: " . ($useFold ? "Fold" : "Superfamily") ." not found for taxon \"$taxon\" in \"$relevanceInfoFileName\";");
	}
	$criterionStr = $taxon2criterionRef->{ $taxon};
	
	# if( $criterionStr =~ m/,.+,/){
	#     print STDERR "ALERT: Using the first ". ($useFold ? "fold" : "superfamily") . " for $taxon from list: $criterionStr!\n";
	# }
	# $criterionStr =~ s/,.*//;
    }
    
    # convert criterion string into a hash
    my %taxonsCriterion;
    foreach my $criterion (split( /,/, $criterionStr)){
	$taxonsCriterion{ $criterion} = 1;
    }
    

    if( $DEBUG){ print STDERR "Using " . join(', ', sort( keys( %taxonsCriterion))) . " for " . ($useFold ? "fold" : "superfamily") . "\n"; }

    $relevantTaxaRef = {};  # initialize
    my $tpsCounter = 0;
    my $dbSize = keys %$taxon2criterionRef;
    
    my $otherTaxon;
    foreach $otherTaxon (sort keys %$taxon2criterionRef){
	my @criteria = split( /,/, $taxon2criterionRef->{ $otherTaxon});
	foreach my $individualCriterion (@criteria){
	    if( exists( $taxonsCriterion{ $individualCriterion})){
		if( $DEBUG){ print STDERR "Found TP: $otherTaxon (" . ($useFold ? "Fold" : "Superfamily") .") -> $individualCriterion (all: $taxon2criterionRef->{$otherTaxon}) \n"; }
		$relevantTaxaRef->{$otherTaxon} = 1;
		$tpsCounter++;
		last;
	    }
	}
    }

    if( $ignoreSelfHit){
	# check for self-hit
	if( defined( $relevantTaxaRef->{ $taxon})){
	    delete( $relevantTaxaRef->{$taxon});
	    $tpsCounter--;
	}
    }
    
    return ($relevantTaxaRef, $tpsCounter, $dbSize);
    
} # END parseTaxaInfo()



sub fileToString 
# $file_ // Input file
{
    my $file_ = $_[0];
    
    my $msg = "Input file '" . $file_ . "' failed to open.\n" . "Died";

    my $terminator = $/;
    undef $/;
    open INPUT, "<$file_" or die $msg;
    my $str = <INPUT>; # terminator undefined : $str is the whole file.
    close INPUT;

    $/ = $terminator; # Restore for normal behavior later

    return $str;
}
    
