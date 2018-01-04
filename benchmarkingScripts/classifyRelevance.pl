#!/usr/bin/env perl

################################################################################
# File: classifyRelevance.pl
# Output: Produces a sorted retrieval list with the columns FP, TP and E-value.
# Input: BLAST or other retrieval stdout file, taxa, superfamily( & fold) definition file
# NOTE: Assumes that the query is in the relevancy list (e.g., is potentially counted in the total number of TPs)
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
# Copyright 2015-2018 Hyrum D. Carroll

# NOTE: The term "Index" in a variable name identifies the location (index value: 0 .... len - 1) of the start, end, etc. of a string or array
# NOTE: the terms "superfamily" and "fold" can be substituted with "more strigent criterion" and "less strigent criterion"

#
# To-Do
# 

use strict;
use warnings;

#use Bio::SearchIO;
#use Bio::SeqIO;

use SortingShuffling qw(sort2Arrays);  # HDC::SortingShuffling

# use Common qw($FALSE $TRUE $VERBOSE $DEBUG VERBOSE DEBUG MIN MAX);
my $FALSE = 0;
my $TRUE = 1;

my $VERBOSE = $FALSE;
my $DEBUG = $FALSE;

sub VERBOSE{
    if( $VERBOSE){ print STDERR "VERBOSE: @_\n"; }
}

sub DEBUG{
    if( $DEBUG){ print STDERR "DEBUG: @_\n"; }
}

sub MIN{
    return ( ($_[0] <= $_[1]) ? $_[0] : $_[1]);
}

sub MAX{
    return ( ($_[0] >= $_[1]) ? $_[0] : $_[1]);
}

my $GAP = '-';

# For the domain information 2D array indicies (for domain locations):
my $DOMAIN_NAME_INDEX = 0;
my $START_INDEX = $DOMAIN_NAME_INDEX + 1;
my $LAST_INDEX = $START_INDEX + 1;

my $QUERY_INDEX = 0;
my $SUBJECT_INDEX = $QUERY_INDEX + 1;
my @TYPE = ('query', 'subject');



# types of hits output files:
my $UNKNOWN_HITS_OUTPUT_TYPE =   -1;
my $BLAST_HITS_OUTPUT_TYPE =     $UNKNOWN_HITS_OUTPUT_TYPE + 1;
my $SIMPLE_HITS_OUTPUT_TYPE =    $BLAST_HITS_OUTPUT_TYPE + 1;
my $GLOBAL_HITS_OUTPUT_TYPE =    $SIMPLE_HITS_OUTPUT_TYPE + 1;
my $PHMMER_HITS_OUTPUT_TYPE =    $GLOBAL_HITS_OUTPUT_TYPE + 1;
my $JACKHMMER_HITS_OUTPUT_TYPE = $PHMMER_HITS_OUTPUT_TYPE + 1;
my $PSISEMIGLOBAL_HITS_OUTPUT_TYPE = $JACKHMMER_HITS_OUTPUT_TYPE + 1;

my $DEFAULT_MIN_OVERLAP = 50;  # percentage (i.e., 50 = 50%)

my $usage = "Usage: $0:
  [--taxon=<taxon name>] (e.g., d12asa_)
  [--family=<superfamily name>] (e.g., 46458 or cd00001)
  [--fold=<fold name>] (e.g., 46457)
  --blastp=<BLAST hits output file> 
  --psiblast=<PSI-BLAST hits output file>
  --psisemiglobal=<PSI-SemiGLOBAL hits output file>
  --simpleHits=<Simplfied BLAST-like hits output file>
  --global=<GLOBAL hits output file>
  --hmmer=<HMMER hits table output file>
  --phmmer=<HMMER hits table output file>
  --jackhmmer=<HMMER hits table output file>
  --rel=<relevance info file (e.g., id2superfamily.tab or id2superfamily_fold.tab)>
  [--spouge=<retrieval list output file>]
  [--spougee=<retrieval list (with --spougeEValue threshold) output file>]
  [--spougeEValue=<threshold E-value for --spougee>]
  [--spougeExt (verbose output in the .spouge file)] 
  [--tab=<retrieval list output file>]
  [--tabe=<retrieval list (with --tabEValue threshold) output file>]
  [--tabEValue=<threshold E-value for --tabe>]
  [--usefold] (use fold (secondary) classification information in --rel for relevance; Default: unset) 
  [--nofold] (Ignore fold (secondary) classification information in --rel; Default: set)
  [--ignoreSelfHit] (Ignore database hits matching --taxon; Default: set if --taxon set)
  [--allowSelfHit] (Allow for all possible database hits; Default: unset unless --taxon set)
  [--randomsAsIrrelevants] (Only sequences in the same superfamily/fold are relevant, taxa labels starting with \"random\" are irrelevants, the rest are ambiguous and therefore ignored; Default: true) 
  [--norandomsAsIrrelevants] (Only sequences in the same superfamily/fold are relevant, the rest are irrelevants; Default: not set)
  [--overlap[=<percentage overlap (e.g., 75)>]] (Default: not used; Default percentage: $DEFAULT_MIN_OVERLAP.  Note: Requires --domainLocs)
  [--domainLocs=<domain locations filename (e.g., family_members-nums-domainLocsOnly.annot)>] (Default: not used)
  [--combineHSPs] (Aggregate all of the HSPs together (for BLAST results)) (Default: not set)
  [--ignoreMultipleHSPs] (Ignore multiple HSPs; Default: set, unless enviroment variable ALLOW_MULTIPLE_HSPS is set) (NOTE: Only implemented for BLAST output)
  [--allowMultipleHSPs] (Allow multiple HSPs (and evaluate each HSP separatly); Default: unset, unless enviroment variable ALLOW_MULTIPLE_HSPS is set)  (NOTE: Only implemented for BLAST output)
  [--ignoreMultipleSeqs] (Ignore subsequent sequences with the same taxon label; Default: set, unless enviroment variable ALLOW_MULTIPLE_SEQS is set)
  [--allowMultipleSeqs] (Allow multiple sequences with the same taxon label; Default: unset, unless enviroment variable ALLOW_MULTIPLE_SEQS is set)
  [--totalTPs=<Number of total relevant records> (Default: calculated from the relevance info file)]
  [-v] (verbose; Default: unset)
  [-d] (debug; Default: unset)";
#  [-out=<output file>]
#  -total=<total number of sequences in the DB>


my %args;  # keys are: 

for( my $argvIndex = 0; $argvIndex <= $#ARGV; $argvIndex++){
    if( $ARGV[ $argvIndex] =~ m/^-?-([^\s=]+)=(.*)/){
	$args{lc($1)} = $2; # lc = lowercase
	#print "$0 found arg $1 (=$args{$1})\n";
	splice @ARGV, $argvIndex, 1;
	$argvIndex--;
    }elsif( $ARGV[ $argvIndex] =~ m/^-?-([^\s]+)/){
	$args{lc($1)} = undef; #$FALSE; # lc = lowercase
	#print "$0 found arg $1 (=$args{$1})\n";
	splice @ARGV, $argvIndex, 1;
	$argvIndex--;
    }
}

my $taxon;
my $superfamily;
my $fold;
my $hitsOutputFileName;
my $hitsOutputFileType = $UNKNOWN_HITS_OUTPUT_TYPE;
my $relevanceInfoFileName;
#my $outputFileName = "";
my $spougeOutputFileName = "";
my $spougeeOutputFileName = "";
my $spougeEValue = "";
my $spougeExtendedOutput = $FALSE;
my $tabOutputFileName = "";
my $tabeOutputFileName = "";
my $tabEValue = "";
#my $totalSeqs = -1;
my $useFold = $FALSE;
my $ignoreSelfHit = $FALSE;   # defaults to $TRUE later if $taxon is defined
my $randomsAsIrrelevants = $TRUE;
my $minimumOverlap = 0;  # 0 means not used
my $domainLocsFilename = "";
my $allowMultipleHSPs = $FALSE;
my $allowMultipleSeqs = $FALSE;  # set for bootstrapping
my $combineHSPs = $FALSE;
my $totalTPs = -1; # NOT USED

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
    print STDERR "\nERROR: Either --taxon, --family or --fold is required!\n\n";
    die( "$usage\n");
}


if( defined($args{"blastp"})){
    $hitsOutputFileName = $args{"blastp"};
    $hitsOutputFileType = $BLAST_HITS_OUTPUT_TYPE;
}elsif( defined($args{"psiblast"})){
    $hitsOutputFileName = $args{"psiblast"};
    $hitsOutputFileType = $BLAST_HITS_OUTPUT_TYPE;
}elsif( defined($args{"psisemiglobal"})){
    if( ! defined( $args{"taxon"})){
	print STDERR "\n\nERROR: --psisemiglobal requires --taxon (due to bug in PSI-SemiGLOBAL output near \"Query=\"\n\n";
	die( "$usage\n");
    }
    $hitsOutputFileName = $args{"psisemiglobal"};
    $hitsOutputFileType = $PSISEMIGLOBAL_HITS_OUTPUT_TYPE;
}elsif( defined($args{"psiblast_fdr"})){
    $hitsOutputFileName = $args{"psiblast_fdr"};
    $hitsOutputFileType = $SIMPLE_HITS_OUTPUT_TYPE;
}elsif( defined($args{"simplehits"})){
    $hitsOutputFileName = $args{"simplehits"};
    $hitsOutputFileType = $SIMPLE_HITS_OUTPUT_TYPE;
}elsif( defined($args{"global"})){
    $hitsOutputFileName = $args{"global"};
    $hitsOutputFileType = $GLOBAL_HITS_OUTPUT_TYPE;
}elsif( defined($args{"hmmer"})){
    $hitsOutputFileName = $args{"hmmer"};
    $hitsOutputFileType = $JACKHMMER_HITS_OUTPUT_TYPE; # maybe should be phmmer?
}elsif( defined($args{"phmmer"})){
    $hitsOutputFileName = $args{"phmmer"};
    $hitsOutputFileType = $PHMMER_HITS_OUTPUT_TYPE;
}elsif( defined($args{"jackhmmer"})){
    $hitsOutputFileName = $args{"jackhmmer"};
    $hitsOutputFileType = $JACKHMMER_HITS_OUTPUT_TYPE;
}elsif( defined($args{"hmmer-non_iterative"})){
    $hitsOutputFileName = $args{"hmmer-non_iterative"};
    $hitsOutputFileType = $PHMMER_HITS_OUTPUT_TYPE;
}elsif( defined($args{"hmmer-iterative"})){
    $hitsOutputFileName = $args{"hmmer-iterative"};
    $hitsOutputFileType = $JACKHMMER_HITS_OUTPUT_TYPE;
}

if( ! defined( $hitsOutputFileName)){
    print STDERR "\n\nERROR: --blastp, --psiblast, --psisemiglobal, --global or --{p|jack}hmmer is required!\n\n";
    die( "$usage\n");
}

if( ! -s $hitsOutputFileName){
    print STDERR "\n\nERROR: hitsOutputFileName \"$hitsOutputFileName\" $!\n";
    print STDERR "pwd: " . `pwd`."\n\n";
    die( "$usage\n");
}


if( defined( $args{"rel"})){
    $relevanceInfoFileName = $args{"rel"};
}
if( ! defined( $relevanceInfoFileName) ||
    $relevanceInfoFileName eq "" ||
    ! -s $relevanceInfoFileName ){
    print STDERR "\n\nERROR: Bad family,superfamily[,fold] info file \"$relevanceInfoFileName\"!\n\n";
    die( "$usage\n");
}

#if( defined( $args{"out"})){
#    $outputFileName = $args{"out"};
#}

if( defined( $args{"spouge"})){
    $spougeOutputFileName = $args{"spouge"};
}else{
    # use STDOUT
}
#if( $spougeOutputFileName eq ""){
#    print STDERR "ERROR: SPOUGE output file not specified!\n";
#    die( "$usage\n");
#}

if( defined( $args{"spougee"})){
    if( ! defined( $args{"spougeevalue"})){
	print STDERR "\nALERT: Ignoring --spougee value since --spougeEValue is not specified!\n\n";
    }else{
	$spougeeOutputFileName = $args{"spougee"};
    }
}

if( defined( $args{"spougeevalue"})){
    if( ! defined( $args{"spougee"})){
	print STDERR "\nALERT: Ignoring --spougeEValue since --spougee is not specified!\n\n";
    }else{
	$spougeEValue = $args{"spougeevalue"};
    }
}

if( exists( $args{"spougeext"}) || $DEBUG){
    $spougeExtendedOutput = $TRUE;
}


if( defined( $args{"tab"})){
    $tabOutputFileName = $args{"tab"};
}else{
    # use STDOUT
}
#if( $tabOutputFileName eq ""){
#    print STDERR "ERROR: TAB output file not specified!\n";
#    die( "$usage\n");
#}

if( defined( $args{"tabe"})){
    if( ! defined( $args{"tabevalue"})){
	print STDERR "\nALERT: Ignoring --tabe value since --tabEValue is not specified!\n\n";
    }else{
	$tabeOutputFileName = $args{"tabe"};
    }
}

if( defined( $args{"tabevalue"})){
    if( ! defined( $args{"tabe"})){
	print STDERR "\nALERT: Ignoring --tabEValue since --tabe is not specified!\n\n";
    }else{
	$tabEValue = $args{"tabevalue"};
    }
}


#if( defined( $args{"tot"})){
#    $totalSeqs = $args{"tot"};
#}elsif( defined( $args{"total"})){
#    $totalSeqs = $args{"total"};
#}
#
#if( $totalSeqs !~ m/^\d+$/ ||
#    $totalSeqs <= 0){
#    print STDERR "ERROR: Non-numerical count of the number of sequences in the database \"$totalSeqs\"!\n";
#    die( "$usage\n");
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

if( exists( $args{"allowselfhit"})){
    $ignoreSelfHit = $FALSE;
}elsif( exists( $args{"ignoreselfhit"})){
    $ignoreSelfHit = $TRUE;
    if(	! defined( $taxon)){
	print STDERR "\n\nERROR: --taxon must be specified to ignore self hits\n\n";
	die( "$usage\n");
    }
}

if( exists( $args{"randomsasirrelevants"})){
    $randomsAsIrrelevants = $TRUE;
}
if( exists( $args{"norandomsasirrelevants"})){
    $randomsAsIrrelevants = $FALSE;
}

if( exists( $args{"overlap"})){
    if( ! defined( $args{"domainlocs"})){
	print STDERR "\n\nERROR: No domain locations file specified! (use --domainLocs=...)\n\n";
	die( "$usage\n");
    }
    
    if( defined( $args{"overlap"})){
	$minimumOverlap = $args{"overlap"};
    }else{
	$minimumOverlap = $DEFAULT_MIN_OVERLAP;
    }
    if( $DEBUG){ print STDERR "HyDEBUG: minimumOverlap: $minimumOverlap\n"; }
}

if( defined( $args{"domainlocs"})){
    $domainLocsFilename = $args{"domainlocs"};
    if( $DEBUG){ print STDERR "HyDEBUG: domainLocsFilename: $domainLocsFilename\n"; }
}


if( exists( $args{"combinehsps"})){
    if( $hitsOutputFileType != $BLAST_HITS_OUTPUT_TYPE &&
	$hitsOutputFileType != $SIMPLE_HITS_OUTPUT_TYPE){
	print STDERR "NOTE: Ignoring --combineHSPs based on the hits output type ($hitsOutputFileType)\n";
	$combineHSPs = $FALSE;
    }else{
	$combineHSPs = $TRUE;
	DEBUG( "Set combineHSPs to $combineHSPs based on --combineHSPs flag\n");
    }
}


if( defined( $ENV{'ALLOW_MULTIPLE_HSPS'})){
    $allowMultipleHSPs = $ENV{'ALLOW_MULTIPLE_HSPS'};
    DEBUG( "Set allowMultipleHSPs to $allowMultipleHSPs based on the enviroment variable ALLOW_MULTIPLE_HSPS\n");
}

if( exists( $args{"allowmultiplehsps"})){
    $allowMultipleHSPs = $TRUE;
    DEBUG( "Set allowMultipleHSPs to $allowMultipleHSPs based on --allowMultipleHSPs flag\n");
}
if( exists( $args{"ignoremultiplehsps"})){
    $allowMultipleHSPs = $FALSE;
    DEBUG( "Set allowMultipleHSPs to $allowMultipleHSPs based on --ignoreMultipleHSPs flag\n");
}


if( defined( $ENV{'ALLOW_MULTIPLE_SEQS'})){
    $allowMultipleSeqs = $ENV{'ALLOW_MULTIPLE_SEQS'};
    DEBUG( "Set allowMultipleSeqs to $allowMultipleSeqs based on the enviroment variable ALLOW_MULTIPLE_SEQS\n");
}

if( exists( $args{"allowmultipleseqs"})){
    $allowMultipleSeqs = $TRUE;
    DEBUG( "Set allowMultipleSeqs to $allowMultipleSeqs based on --allowMultipleSeqs flag\n");
}
if( exists( $args{"ignoremultipleseqs"})){
    $allowMultipleSeqs = $FALSE;
    DEBUG( "Set allowMultipleSeqs to $allowMultipleSeqs based on --ignoreMultipleSeqs flag\n");
}

if( defined( $args{'totaltps'})){
    $totalTPs = $args{'totaltps'};
    DEBUG( "Using $totalTPs as the total number of relevant records\n");
}


if( $DEBUG){
    print STDERR "$0:\n";
    print STDERR "\ttaxon: $taxon\n";
    if( defined( $superfamily)){
	print STDERR "\tsuperfamily: $superfamily\n";
    }
    if( defined( $fold)){
	print STDERR "\tfold: $fold\n";
    }
    print STDERR "\thitsOutputFileName: $hitsOutputFileName\n";
    print STDERR "\thitsOutputFileType: $hitsOutputFileType\n";
    print STDERR "\trelevanceInfoFileName: $relevanceInfoFileName\n";
    print STDERR "\tspougeOutputFileName: $spougeOutputFileName\n";
    print STDERR "\tspougeeOutputFileName: $spougeeOutputFileName\n";
    print STDERR "\tspougeEValue: $spougeEValue\n";
    print STDERR "\tspougeExtendedOutput: $spougeExtendedOutput\n";
    print STDERR "\ttabOutputFileName: $tabOutputFileName\n";
    print STDERR "\ttabeOutputFileName: $tabeOutputFileName\n";
    print STDERR "\ttabEValue: $tabEValue\n";
    print STDERR "\tuseFold: $useFold\n";
    print STDERR "\tignoreSelfHit: $ignoreSelfHit\n";
    print STDERR "\trandomsAsIrrelevants: $randomsAsIrrelevants\n";
    print STDERR "\tminimumOverlap: $minimumOverlap\n";
    print STDERR "\tdomainLocsFilename: $domainLocsFilename\n";
    print STDERR "\tallowMultipleHSPs: $allowMultipleHSPs\n";
    print STDERR "\tallowMultipleSeqs: $allowMultipleSeqs\n";
    print STDERR "\tcombineHSPs: $combineHSPs\n";
    print STDERR "\ttotalTPs: $totalTPs\n";
    print STDERR "\tVERBOSE: $VERBOSE\n";
    print STDERR "\tDEBUG: $DEBUG\n";
}


#
# Read in relevant family members
#

my ($relevantTaxaRef, $tpsCounter, $totalSeqs) = parseTaxaInfo($relevanceInfoFileName, $taxon, $criterionStr, $ignoreSelfHit);
DEBUG( "HyDEBUG: tpsCounter: $tpsCounter (ignoreSelfHit: $ignoreSelfHit)");
if( $totalTPs >= 0){
    DEBUG( "Ignoring the total number of relevant records from the relevance info file and using the valued specified on the command-line, $totalTPs");
    $tpsCounter = $totalTPs;
}
my $totalNumFPs = $totalSeqs - $tpsCounter + $ignoreSelfHit;

#
# parse and store domain locations
#
my %domainIndices;  # keys are taxon labels; values are 2D arrays with the domain name, start, and end as columns
my $calculateOverlaps = $FALSE;
if( $domainLocsFilename ne "" &&
    -s $domainLocsFilename &&
    $minimumOverlap > 0){

    parseDomainLocs( $domainLocsFilename, \%domainIndices);
    if( ! exists( $domainIndices{ $taxon})){
	print STDERR "ERROR: $taxon not found in the domainIndices hash!\n";
	print STDERR "%domainIndices: %domainIndices\n";
	exit(1);
    }

    $calculateOverlaps = $TRUE;
}


# keep track of the E-value to calculate the ROC
my @evalues;
my @hitLabels;


my $TAXON_LABEL            = 0;
#my $EVALUE                 = 1;
my $OVERLAP                = 1;  # string of aligned residues (only if calculating the overlap) 


my $data = [[]]; # Index 0: hit index; index 1: TAXON_LABEL, 
my @eValues;
my @origIndices;

if( $hitsOutputFileType == $BLAST_HITS_OUTPUT_TYPE){
    parseBlastHits($data, \@eValues, \@origIndices, $hitsOutputFileName, $ignoreSelfHit, $calculateOverlaps, $allowMultipleSeqs, $allowMultipleHSPs, $combineHSPs);
}elsif( $hitsOutputFileType == $SIMPLE_HITS_OUTPUT_TYPE){
    parseSimpleHits($data, \@eValues, \@origIndices, $hitsOutputFileName, $calculateOverlaps, $allowMultipleSeqs);
}elsif( $hitsOutputFileType == $GLOBAL_HITS_OUTPUT_TYPE){
    parseGlobalHits($data, \@eValues, \@origIndices, $hitsOutputFileName, $totalSeqs, $calculateOverlaps, $allowMultipleSeqs);
}elsif( $hitsOutputFileType == $PHMMER_HITS_OUTPUT_TYPE || $hitsOutputFileType == $JACKHMMER_HITS_OUTPUT_TYPE){
    parseHmmerHits($data, \@eValues, \@origIndices, $hitsOutputFileName, $ignoreSelfHit, $calculateOverlaps, $allowMultipleSeqs, $allowMultipleHSPs, $combineHSPs, $taxon);
    #parseHmmerHits($data, \@eValues, \@origIndices, $hitsOutputFileName, $totalSeqs, $calculateOverlaps, $allowMultipleSeqs);
}elsif( $hitsOutputFileType == $PSISEMIGLOBAL_HITS_OUTPUT_TYPE){
    parsePsiSemiGlobalHits($data, \@eValues, \@origIndices, $hitsOutputFileName, $ignoreSelfHit, $calculateOverlaps, $taxon, $allowMultipleSeqs);
}else{
    die( "\n\nERROR: Unknown hits output type: $hitsOutputFileType!");
}

#
# Output tab and ROC files
#

my %reportedHitLabels;

#	if( $outputFileName ne ""){
#	    $outputFileName    .= ".$iter";
#	}
#}
#
#if( $outputFileName ne ""){
#	open(OUT, ">$outputFileName") or die( "ERROR: Can not open / create output file \"$outputFileName\": $! ");
#	#print OUT "#Pfam_family\tprobe\tprobe_start_index\thit_accession\thit_start_index\tE-value\thit_type(s)\toverlap_percentage(s)\tcoverage_percentage(s)\n";
#	print OUT "#taxon\tprobe\thit_accession\tE-value\n";
my $spougeHeader = "";
if( $spougeOutputFileName ne ""){
    $spougeHeader = $spougeOutputFileName;
}elsif( defined( $taxon)){
    $spougeHeader = $taxon;
}elsif( defined( $superfamily)){
    $spougeHeader = $superfamily;
}elsif( defined( $fold)){
    $spougeHeader = $fold;
}
$spougeHeader .= "\n";
$spougeHeader .= "$tpsCounter\n";

if( $spougeOutputFileName ne ""){
    open(SPOUGE, ">$spougeOutputFileName") or die( "ERROR: Can not open / create SPOUGE output file \"$spougeOutputFileName\": $! ");
    print SPOUGE "$spougeHeader";
}

if( $spougeeOutputFileName ne ""){
    open(SPOUGEE, ">$spougeeOutputFileName") or die( "ERROR: Can not open / create SPOUGEe output file \"$spougeeOutputFileName\": $! ");
    print SPOUGEE "$spougeHeader";
}


my $tabHeader = "#FP($totalNumFPs)\tTP($tpsCounter)\tE-value\tFPR($totalNumFPs)\tTPR($tpsCounter)\n";
if( $tabOutputFileName ne ""){
    open(TAB, ">$tabOutputFileName") or die( "ERROR: Can not open / create TAB output file \"$tabOutputFileName\": $! ");
    print TAB "$tabHeader";
}

if( $tabeOutputFileName ne ""){
    open(TABE, ">$tabeOutputFileName") or die( "ERROR: Can not open / create TABe output file \"$tabeOutputFileName\": $! ");
    print TABE "$tabHeader";
}


my $tabTpCounter = 0;
my $tabFpCounter = 0;

# NOTE: tabe just uses tab[TF]pCounter

for( my $i = 0; $i <= $#origIndices; $i++){

    my $j = $origIndices[$i];
    my $hitLabel           = $data->[$j][$TAXON_LABEL];
    my $eValue             = $eValues[$i];  # use "$i" b/c @eValues is already sorted
    my $relevantRecord     = $FALSE;
    my $tabTP              = $FALSE;

    if( $DEBUG){ print STDERR "DEBUG:\t$hitLabel:"; }

    if( $DEBUG){
	if( exists( $reportedHitLabels{$hitLabel})){
	    DEBUG( "\tFound duplicate taxon label ($reportedHitLabels{$hitLabel} times)");
	}
	$reportedHitLabels{$hitLabel}++;
    }

    if( exists($relevantTaxaRef->{$hitLabel})){
	# relevant record
	$relevantRecord = $TRUE;
	if( $DEBUG){ print STDERR " names match\n"; }
	
	if( $calculateOverlaps){
	    # determine if overlap meets minimum threshold
	    if( $DEBUG){ print STDERR "OVERLAP_DEBUGGING:\tdata->[$j][$OVERLAP]: $data->[$j][$OVERLAP]\n"; }
	    if( $data->[$j][$OVERLAP] >= $minimumOverlap){
		$tabTP = $TRUE;
	    }
	}else{
	    $tabTP = $TRUE;
	}

	#if( $outputFileName ne ""){
	#    print OUT "$taxon\t$probe\t$hitLabel\t$eValue\n";
	#}
    }elsif( $ignoreSelfHit == $TRUE &&
	    $hitLabel eq $taxon ){
	if( $DEBUG){ print STDERR "DEBUG:\t\tIgnoring self hit: $hitLabel\n"; }
	next;
    }elsif( $TRUE && $DEBUG){
	# irrelevant record
	print STDERR " not part of ".($useFold ? "Fold" : "Superfamily") ."\n";
    }

    if( $tabTP){
	$tabTpCounter++;
    }elsif( $randomsAsIrrelevants && $hitLabel !~ m/^RANDOM/){
	# not part of the same superfamily/fold and not a "random\d{6}" sequence, so ignoring 
	next;
    }else{
	$tabFpCounter++;
    }
    
    my $fpr = 0;
    my $tpr = 0;
    if ($totalNumFPs != 0){
	$fpr=$tabFpCounter/$totalNumFPs;
    }
    if( $tpsCounter != 0){
	$tpr = $tabTpCounter/$tpsCounter;
    }

    my $spougeStr = "$tabTP\t$eValue";
    if( $spougeExtendedOutput == $TRUE){
	$spougeStr .= "\t\t$hitLabel"; # leave a blank spot for the weight
	if( $calculateOverlaps){
	    $spougeStr .= "\t$data->[$j][$OVERLAP]";
	}
	$spougeStr .= "\t$relevantRecord";
    }
    $spougeStr .= "\n";
    
    if( $spougeOutputFileName ne ""){
	print SPOUGE $spougeStr;
    }

    if( $spougeeOutputFileName ne "" &&
	$eValue <= $spougeEValue){
	print SPOUGE $spougeStr;
    }
    
    my $tabStr = sprintf "$tabFpCounter\t$tabTpCounter\t$eValue\t%.7f\t%.7f\n", $fpr, $tpr;
    if( $tabOutputFileName ne ""){
	print TAB "$tabStr";
    }

    if( $tabeOutputFileName ne "" &&
	$eValue <= $tabEValue){
	print TABE "$tabStr";
    }
} # END hit loop


#if( $outputFileName ne ""){
#	close( OUT);
#}

if( $spougeOutputFileName ne ""){
    close( SPOUGE);
}
if( $spougeeOutputFileName ne ""){
    close( SPOUGEE);
}

if( $tabOutputFileName ne ""){
    close( TAB);
}
if( $tabeOutputFileName ne ""){
    close( TABE);
}

# if( $foundMultipleHitsPerSequence == $TRUE &&
#     $NEVER_DO_FULL_FILES == $FALSE){
#     print STDERR "\nALERT: Multiple hits per database sequence not supported!\n\n";
# }

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
		#    print STDERR "\ALERT: Ignoring previous fold ($taxon2criterionRef->{uc($1)}) for $1 in favor of $3!\n\n";
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
		#    print STDERR "\nALERT: Ignoring previous superfamily ($taxon2criterionRef->{uc($1)}) for $1 in favor of $2!\n\n";
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
	#     print STDERR "\nALERT: Using the first ". ($useFold ? "fold" : "superfamily") . " for $taxon from list: $criterionStr!\n\n";
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
	    if( $DEBUG){ print STDERR "Removing self-hit from relevant records list\n"; }
	}
    }
    
    return ($relevantTaxaRef, $tpsCounter, $dbSize);
    
} # END parseTaxaInfo()



#                  alignmentSeq starts here
#                  |              residueIndexToFind example1
#                  |              |             residueIndexToFind example2
#                  |              |             |
#                  v              v             v
# ...rrrrrrrrrrrrrrrrrrrrrrrr---rrrrrrrrrrrr--rrrrrrrr      (where "r" is any residue and "-" is a gap)
# sequence index:  0              12            24   29   # indices relative to the residues in the string
# alignment index: 0         10   15            29   34   # indices relative to the alignment string
# residue index:   100            112           124  129  # indices relative to the query/subject sequence (not used by getAlignmentIndex())
#
# example1: getAlignmentIndex( alignmentSeq: "rrrrrrrrrr---rrrrrrrrrrrr--rrrrrrrr",
#                              residueIndexToFind: 12, 
#                              residueIndex: 0 (in); 12 (out), 
#                              alignmentIndex: 0 (in); 15 (out))
#
# example2: getAlignmentIndex( alignmentSeq: "rrrrrrrrrr---rrrrrrrrrrrr--rrrrrrrr",
#                              residueIndexToFind: 24, 
#                              residueIndex: 12 (in); 24 (out), 
#                              alignmentIndex: 15 (in); 29 (out))

# Essentially calculates the number of gaps found in $seq between  
sub getAlignmentIndex{
    # alignmentSeq: in
    # residueIndexToFind: in
    # residueIndex: in/out
    # alignmentIndex: in/out
    # outOfBounds: out  ('<' if residueIndexToFind is before residueIndexRef;
    #                    '>' if residueIndexToFind is after (greater than) the length of alignmentSeq;
    #                    if set, residueIndex and alignmentIndex are ... )
    my ( $alignmentSeq, $residueIndexToFind, $residueIndexRef, $alignmentIndexRef, $outOfBoundsRef) = @_;
    DEBUG( "getAlignmentIndex( $alignmentSeq, $residueIndexToFind, $$residueIndexRef, $$alignmentIndexRef, $$outOfBoundsRef)");

    # check for out of bounds cases:
    if( $residueIndexToFind < $$residueIndexRef){
	DEBUG("\nALERT: residueIndexToFind ($residueIndexToFind) < residueIndexRef ($$residueIndexRef)!\n");
	$$outOfBoundsRef = '<';
        # $$residueIndexRef  # no change
	# $$alignmentIndexRef  # no change
	return;
    }
    
    my $alignmentSeqLen = length( $alignmentSeq);
    if( $residueIndexToFind >= $alignmentSeqLen){
	DEBUG("\nALERT: residueIndexToFind ($residueIndexToFind) > alignmentSeqLen ($alignmentSeqLen)!\n");
	$$outOfBoundsRef = '>';
	# return the largest index in the sequence
	$$alignmentIndexRef = $alignmentSeqLen - 1;
	return;
    }

    # # if the residue index we're looking for is actually before the starting point (residueIndexRef), than it will be assigned to the starting point
    # if( $$residueIndexRef != 0  &&  $residueIndexToFind < $$residueIndexRef){ print STDERR "\nALERT: residueIndexToFind: $residueIndexToFind) < residueIndexRef: $$residueIndexRef)!\n\n"; }
    
    my @alignmentSeq = split(//, $alignmentSeq);
    # start with the residueIndexRef, and increment alignmentIndex when there's NOT a gap
    while( $$residueIndexRef < $residueIndexToFind  && $$alignmentIndexRef < $alignmentSeqLen ){
	if( $alignmentSeq[ $$alignmentIndexRef ] ne $GAP){
	    $$residueIndexRef++;
	}
	$$alignmentIndexRef++;
    }
    DEBUG( "\t residueIndex: $$residueIndexRef, alignmentIndexRef: $$alignmentIndexRef");
} # END getAlignmentIndex()


# Get the intersection of the query and subject domains
# Returns the number of residues in the query domain(s)
sub getQuerySubjectDomainIntersection{
    my ($domainIndicesRef, $queryLabel, $subjectLabel, $domainStartIndices, $domainLastIndices) = @_;

    my @labels = ($queryLabel, $subjectLabel);
    
    my @domainsRefs;
    for( my $typeI = 0; $typeI <= $#TYPE; $typeI++){
	if( ! defined( $domainIndicesRef->{ $labels[$typeI]})){
	    die( "ERROR: Domain location not found for $labels[$typeI]!");
	}
	push( @domainsRefs, $domainIndicesRef->{ $labels[$typeI] });
	if( $DEBUG){
	    my $str = "$labels[$typeI] domains (".($#{$domainsRefs[$typeI]} + 1)."):";
	    for( my $domainInfoI = 0; $domainInfoI <= $#{$domainsRefs[$typeI]}; $domainInfoI++){
		$str .= " $domainsRefs[$typeI]->[$domainInfoI][$DOMAIN_NAME_INDEX]";
	    }
	    DEBUG( "\t" . $str);
	}
    }
    my $numQueryDomainResidues = 0;
    for( my $domainInfoQueryI = 0; $domainInfoQueryI <= $#{$domainsRefs[$QUERY_INDEX]}; $domainInfoQueryI++){
	my $domainName = $domainsRefs[$QUERY_INDEX]->[$domainInfoQueryI][$DOMAIN_NAME_INDEX];
	# look for domainName for SUBJECT_INDEX
	for( my $domainInfoSubjectI = 0; $domainInfoSubjectI <= $#{$domainsRefs[$SUBJECT_INDEX]}; $domainInfoSubjectI++){
	    if( $domainsRefs[$SUBJECT_INDEX]->[$domainInfoSubjectI][$DOMAIN_NAME_INDEX] eq $domainName){
		push( @{$domainStartIndices->[$QUERY_INDEX]}, $domainsRefs[$QUERY_INDEX]->[$domainInfoQueryI][$START_INDEX]);
		push( @{$domainLastIndices->[$QUERY_INDEX]},  $domainsRefs[$QUERY_INDEX]->[$domainInfoQueryI][$LAST_INDEX]);
		
		push( @{$domainStartIndices->[$SUBJECT_INDEX]}, $domainsRefs[$SUBJECT_INDEX]->[$domainInfoSubjectI][$START_INDEX]);
		push( @{$domainLastIndices->[$SUBJECT_INDEX]},  $domainsRefs[$SUBJECT_INDEX]->[$domainInfoSubjectI][$LAST_INDEX]);
		
		last;
	    }
	}
	$numQueryDomainResidues += $domainsRefs[$QUERY_INDEX]->[$domainInfoQueryI][$LAST_INDEX] - $domainsRefs[$QUERY_INDEX]->[$domainInfoQueryI][$START_INDEX] + 1;
    }
    DEBUG( "\tnumQueryDomainResidues: $numQueryDomainResidues");
    if( $#{$domainStartIndices->[$QUERY_INDEX]} != $#{$domainStartIndices->[$SUBJECT_INDEX]}){
	die( "$#{$domainStartIndices->[$QUERY_INDEX]} != $#{$domainStartIndices->[$SUBJECT_INDEX]}");
    }
    if( $#{$domainLastIndices->[$QUERY_INDEX]} != $#{$domainLastIndices->[$SUBJECT_INDEX]}){
	die( "$#{$domainLastIndices->[$QUERY_INDEX]} != $#{$domainLastIndices->[$SUBJECT_INDEX]}");
    }

    return $numQueryDomainResidues;
}


sub printDomains{
    my ( $seqsRef, $alignmentStartIndex, $alignmentLastIndex, $indent) = @_;
    #DEBUG("printDomains( ($seqsRef->[$QUERY_INDEX] (".length($seqsRef->[$QUERY_INDEX])."), $seqsRef->[$SUBJECT_INDEX] (".length($seqsRef->[$SUBJECT_INDEX]).")), $alignmentStartIndex, $alignmentLastIndex)");
    if( ! defined( $indent)){
	$indent = '';
    }
    DEBUG( "${indent}Printing out domains (as they appear in the alignment) . . .");
    DEBUG( $indent . substr( $seqsRef->[$QUERY_INDEX], $alignmentStartIndex, $alignmentLastIndex - $alignmentStartIndex + 1));
    DEBUG( $indent . substr( $seqsRef->[$SUBJECT_INDEX], $alignmentStartIndex, $alignmentLastIndex - $alignmentStartIndex + 1));
}

sub printDomains4{
    my ( $seqsRef, $queryResidueStartIndex, $queryResidueLastIndex, $subjectResidueStartIndex, $subjectResidueLastIndex) = @_;
    DEBUG("printDomains( ($seqsRef->[$QUERY_INDEX], $seqsRef->[$SUBJECT_INDEX]), $queryResidueStartIndex, $queryResidueLastIndex, $subjectResidueStartIndex, $subjectResidueLastIndex)");
    DEBUG( "Printing out domains (as they appear in the alignment) . . .");
    my $space = ' ' x  MAX( 0, $queryResidueStartIndex - $subjectResidueStartIndex);
    DEBUG( $space . substr( $seqsRef->[$QUERY_INDEX], $queryResidueStartIndex, $queryResidueLastIndex - $queryResidueStartIndex + 1) . " (alignmentIndices: $queryResidueStartIndex..$queryResidueLastIndex)");
    $space = ' ' x  MAX( 0, $subjectResidueStartIndex - $queryResidueStartIndex);
    DEBUG( $space . substr( $seqsRef->[$SUBJECT_INDEX], $subjectResidueStartIndex, $subjectResidueLastIndex - $subjectResidueStartIndex + 1) . " (alignmentIndices: $subjectResidueStartIndex..$subjectResidueLastIndex)");
}



# Calculate the percentage of overlap of the alignment of the domains in the query and the subject sequences
sub calculateMultiHspOverlap{
    my ($taxon, $hitLabel, $hspsRef, $domainIndicesRef, $queryDomainNamesRef) = @_;
    # Foreach HSP
    #   Foreach query domain
    #     Calculate the first and last indices of the query domain, that are part of the alignment (shortening the range is necessary)
    #     Map those query alignment indicies to subject alignment indices (trival, they're the same indicies)
    #     Foreach subject domain
    #       If the subject domain name matches the query domain name
    #         Calculate the first and last indices of the subject domain, that are part of the alignment (shortening the range is necessary)
    #         If the range is not empty
    #           Map the (possibly shortened) subject domain indices back to the query sequence (via the alignment again)
    #           Calculate the overlap for this query and subject domain pair
    #           Save the overlap if it's longer than any that were previously found for the this query and subject domain pair
    # Until all query domains are included, (greedily) find the longest overlap for between any pair of unused query and subject domains
    # Sum the overlaps found in the previous step

    # Example 1:
    #             +------- query domain -------+
    #             |                            |
    #             |      3--adj.overlap-4      |
    #             |      |              |      |
    #             v      v              v      v
    # Query:   XXXXXXX---X..............X---XXXXXXX
    # Subject: XX---XXXXXX..............XXXXXX---XX
    #               ^  ^                  ^  ^
    #               |  |                  |  |
    #               1  +- subject domain -+  2
    #
    # Event 1: The start of the query   domain aligns to a gap on the subject sequence (so cooresponding position for the subject sequence is adjusted (higher)
    # Event 2: The end   of the query   domain aligns to a gap on the subject sequence (so cooresponding position for the subject sequence is adjusted (lower)
    # Event 3: The start of the subject domain aligns to a gap on the query   sequence (so cooresponding position for the subject sequence is adjusted (higher)
    # Event 4: The end   of the subject domain aligns to a gap on the query   sequence (so cooresponding position for the subject sequence is adjusted (lower)
    #
    # Example 2:
    #             +------- query domain --+
    #             |                       |
    #             |      +--adj.overlap---+
    #             |      |                |
    #             v      v                v
    # Query:   XXXXXXX---X..............XXXXXX
    # Subject: XX---XXXXXX..............XXXXXX
    #                  ^                  ^  ^
    #                  |                  |  |
    #                  |                  5  |
    #                  |                     |
    #                  +- subject domain-----+
    #
    # Note 5: Overlap is NOT extended to the last subject domain residue
    

    DEBUG("calculateMultiHspOverlap( $taxon, $hitLabel, hspsRef, domainIndicesRef, queryDomainNamesRef)");
    my @labels = (uc($taxon), uc($hitLabel));

    for( my $typeI = 0; $typeI <= $#TYPE; $typeI++){
	if( ! defined( $domainIndicesRef->{ $labels[$typeI]})){
	    die( "ERROR: Domain location not found for $labels[$typeI]!");
	}
    }

    my $numberOfQueryDomains   = $#{$domainIndicesRef->{ $labels[$QUERY_INDEX] }} + 1;
    my $numberOfSubjectDomains = $#{$domainIndicesRef->{ $labels[$SUBJECT_INDEX] }} + 1;

    #
    # Quick check that there's at least one matching domain name between the query and the subject
    #
    my $subjectDomainIndex = 0;
    for( $subjectDomainIndex = 0; $subjectDomainIndex < $numberOfSubjectDomains; $subjectDomainIndex++){
	my $subjectDomainName = $domainIndicesRef->{ $labels[$SUBJECT_INDEX] }[$subjectDomainIndex][$DOMAIN_NAME_INDEX];
	if( exists( $queryDomainNamesRef->{ $subjectDomainName }) ){
	    last;
	}
    }

    if( $subjectDomainIndex >= $numberOfSubjectDomains){
	DEBUG( "\tNo matching domain names between the query and the subject");
	my $overlap = 0;
	return $overlap;
    }
    
    my $domainOverlaps = [[]];  # 2D array: rows: best overlap value for query domains;  columns: best overlap value for subject domains
    for( my $queryDomainIndex = 0; $queryDomainIndex < $numberOfQueryDomains; $queryDomainIndex++){
	for( my $subjectDomainIndex = 0; $subjectDomainIndex < $numberOfSubjectDomains; $subjectDomainIndex++){
	    $domainOverlaps->[$queryDomainIndex][$subjectDomainIndex] = 0;
	}
    }
    
    # foreach HSP
    for( my $hspI = 0; $hspI <= $#{$hspsRef}; $hspI++){
	my $hsp = $hspsRef->[$hspI];

	calculateHspOverlap( $taxon, $hitLabel, $hsp, $domainIndicesRef, $domainOverlaps);
    } # END foreach HSP
    
    # Greedily choose domain pairs based on the longest overlap
    my $queryOverlappingResidues = 0;  # number of query domain residues in the alignment (where the query and the subject domains match)
    my @queryDomainsUsed   = ($FALSE) x $numberOfQueryDomains;   # query   domains already used
    my @subjectDomainsUsed = ($FALSE) x $numberOfSubjectDomains; # subject domains already used
    my $currMaxOverlap = 0;
    do{
	$currMaxOverlap = 0;
	my $currMaxQueryIndex = -1;
	my $currMaxSubjectIndex = -1;
	for( my $queryDomainIndex = 0; $queryDomainIndex < $numberOfQueryDomains; $queryDomainIndex++){
	    if( $queryDomainsUsed[$queryDomainIndex] == $TRUE){
		next; # skip over query domains already choosen 
	    }
	    for( my $subjectDomainIndex = 0; $subjectDomainIndex < $numberOfSubjectDomains; $subjectDomainIndex++){
		if( $subjectDomainsUsed[$subjectDomainIndex] == $TRUE){
		    next; # skip over subject domains already choosen 
		}
		if( $domainOverlaps->[$queryDomainIndex][$subjectDomainIndex] > $currMaxOverlap){
		    $currMaxOverlap = $domainOverlaps->[$queryDomainIndex][$subjectDomainIndex];
		    $currMaxQueryIndex   = $queryDomainIndex;
		    $currMaxSubjectIndex = $subjectDomainIndex;
		}
	    }# END foreach subject domain
	} # END foreach query domain

	if( $currMaxOverlap > 0){
	    $queryDomainsUsed[  $currMaxQueryIndex]   = $TRUE;
	    $subjectDomainsUsed[$currMaxSubjectIndex] = $TRUE;
	    $queryOverlappingResidues += $currMaxOverlap;
	    DEBUG("\tmax overlap: $currMaxOverlap (query: $currMaxQueryIndex, subject: $currMaxSubjectIndex); queryOverlappingResidues: $queryOverlappingResidues");
	}

    }while( $currMaxOverlap > 0 );

    
    my $queryTotalResidues = 0; # the number of residues in all of the query domains 
    for( my $queryDomainIndex = 0; $queryDomainIndex < $numberOfQueryDomains; $queryDomainIndex++){
	my $queryStartIndex = $domainIndicesRef->{ $labels[$QUERY_INDEX] }[$queryDomainIndex][$START_INDEX];
	my $queryLastIndex  = $domainIndicesRef->{ $labels[$QUERY_INDEX] }[$queryDomainIndex][$LAST_INDEX];
	
	$queryTotalResidues += $queryLastIndex - $queryStartIndex + 1;
    }	    

    my $overlap = ($queryOverlappingResidues / $queryTotalResidues) * 100.0;
    DEBUG( "\toverlap: $overlap\% (($queryOverlappingResidues / $queryTotalResidues) * 100.0)\n");
    if( $overlap < 0 ||
	$overlap > 100){
	die("HSP overlap: $overlap!");
    }
    
    return $overlap;
} # END  calculateMultiHspOverlap()


# Calculate the percentage of overlap of the alignment of the domains in the query and the subject sequences
sub calculateHspOverlap{
    my ($taxon, $hitLabel, $hsp, $domainIndicesRef, $domainOverlapsRef) = @_;

    DEBUG("calculateHspOverlap( $taxon, $hitLabel, hsp (E-value: ".$hsp->evalue."), domainIndicesRef" . (defined( $domainOverlapsRef) ? ", domainOverlaps" : "") . ")");
    my @labels = (uc($taxon), uc($hitLabel));

    for( my $typeI = 0; $typeI <= $#TYPE; $typeI++){
	if( ! defined( $domainIndicesRef->{ $labels[$typeI]})){
	    die( "ERROR: Domain location not found for $labels[$typeI]!");
	}
    }

    my @startIndices;
    my @lastIndices;
    for( my $typeI = 0; $typeI <= $#TYPE; $typeI++){
	$startIndices[$typeI] = $hsp->start($TYPE[$typeI]) - 1; # "- 1" to convert from position to index
	$lastIndices[$typeI]  = $hsp->end($TYPE[$typeI]) - 1; # "- 1" to convert from position to index
	# DEBUG( "\tstartIndices[$TYPE[$typeI]]: $startIndices[$typeI], lastIndices[$TYPE[$typeI]]: $lastIndices[$typeI]");
    }
    
    my @alignmentSeqs; # indices: query, subject
    foreach my $seq ($hsp->get_aln()->each_seq()){
	push( @alignmentSeqs, $seq->seq());
    }

    my $alignmentTotalLength = $hsp->length('total');
    
    return calculateAlignmentOverlap( $taxon, $hitLabel, $alignmentTotalLength, \@alignmentSeqs, \@startIndices, \@lastIndices, $domainIndicesRef, $domainOverlapsRef);
} # END calculateHspOverlap()


sub calculateAlignmentOverlap{
    my ($taxon, $hitLabel, $alignmentTotalLength, $alignmentSeqsRef, $startIndicesRef, $lastIndicesRef, $domainIndicesRef, $domainOverlapsRef) = @_;

    my @labels = (uc($taxon), uc($hitLabel));

    for( my $typeI = 0; $typeI <= $#TYPE; $typeI++){
	DEBUG( "\tstartIndicesRef->[$TYPE[$typeI]]: $startIndicesRef->[$typeI], lastIndicesRef->[$TYPE[$typeI]]: $lastIndicesRef->[$typeI]");
    }
    
    DEBUG("\tQuery (".length($alignmentSeqsRef->[$QUERY_INDEX])."): $alignmentSeqsRef->[$QUERY_INDEX]");
    DEBUG("\tSubje (".length($alignmentSeqsRef->[$SUBJECT_INDEX])."): $alignmentSeqsRef->[$SUBJECT_INDEX]");

    # Assuming that both a gap never aligns with a gap
    
    # index 2 index example 1:   
    #          0         1  Result:
    #          01234567890                               01234567890
    # Query:   XXXXXXX--XX      queryIndex2SubjectIndex: 011112367   
    # Subject: XX---XXXXXX      subjectIndex2QueryIndex: 01566678

    # NOTE: If the starting index is not 0, the starting index is added to every element of the respective array.
    # NOTE: When there's a gap in the other sequence, the array points to the residue at the beginning of the gap.  This naturally accounts for the alignments at the end of the domain.  For the start of a domain (and there's a gap in the other sequence), an adjustment of +1 needs to be made.

    my @queryIndex2SubjectIndex;
    my @subjectIndex2QueryIndex;
    my @queryIndicesAlignsToAGap = (0) x $alignmentTotalLength;    # indices are sequence index - starting index;  values: 1 for a gap, otherwise 0 
    my @subjectIndicesAlignsToAGap = (0) x $alignmentTotalLength;  # indices are sequence index - starting index;  values: 1 for a gap, otherwise 0 
    my $queryGaps = 0;
    my $subjectGaps = 0;
    my @queryAlignment   = split( //, $alignmentSeqsRef->[$QUERY_INDEX]);
    my @subjectAlignment = split( //, $alignmentSeqsRef->[$SUBJECT_INDEX]);

    
    for( my $i = 0; $i <= $#queryAlignment; $i++){
	
	if( $queryAlignment[$i]   eq $GAP){ $queryGaps++;   $subjectIndicesAlignsToAGap[$i - $subjectGaps] = 1; }
	if( $subjectAlignment[$i] eq $GAP){ $subjectGaps++; $queryIndicesAlignsToAGap[  $i - $queryGaps]   = 1; }

	if( $queryAlignment[$i] ne $GAP){
	    $queryIndex2SubjectIndex[$i - $queryGaps]   = $i + $startIndicesRef->[$SUBJECT_INDEX] - $subjectGaps;
	}
	
	if( $subjectAlignment[$i] ne $GAP){
	    $subjectIndex2QueryIndex[$i - $subjectGaps] = $i + $startIndicesRef->[$QUERY_INDEX]   - $queryGaps;
	}
    }
    DEBUG("\t\@queryIndex2SubjectIndex (".($#queryIndex2SubjectIndex + 1)."): @queryIndex2SubjectIndex");
    DEBUG("\t\@subjectIndex2QueryIndex (".($#subjectIndex2QueryIndex + 1)."): @subjectIndex2QueryIndex");

    DEBUG("\t\@queryIndicesAlignsToAGap   (".($#queryIndicesAlignsToAGap + 1)."): @queryIndicesAlignsToAGap");
    DEBUG("\t\@subjectIndicesAlignsToAGap (".($#subjectIndicesAlignsToAGap + 1)."): @subjectIndicesAlignsToAGap");

    my $numberOfQueryDomains   = $#{$domainIndicesRef->{ $labels[$QUERY_INDEX] }} + 1;
    my $numberOfSubjectDomains = $#{$domainIndicesRef->{ $labels[$SUBJECT_INDEX] }} + 1;
    
    my $queryTotalResidues = 0;
    my $queryOverlappingResidues = 0;
    
    for( my $queryDomainIndex = 0; $queryDomainIndex < $numberOfQueryDomains; $queryDomainIndex++){
	my $queryDomainName = $domainIndicesRef->{ $labels[$QUERY_INDEX] }[$queryDomainIndex][$DOMAIN_NAME_INDEX];
	my $queryStartIndex = $domainIndicesRef->{ $labels[$QUERY_INDEX] }[$queryDomainIndex][$START_INDEX];
	my $queryLastIndex  = $domainIndicesRef->{ $labels[$QUERY_INDEX] }[$queryDomainIndex][$LAST_INDEX];
	
	my $queryDomainResidues += $queryLastIndex - $queryStartIndex + 1;
	$queryTotalResidues += $queryDomainResidues;
	DEBUG("\tqueryDomainIndex: $queryDomainIndex ($queryDomainName: $queryStartIndex..$queryLastIndex) (queryDomainResidues: $queryDomainResidues)");

	
	# Adjust the first residue based on the starting position and the start of the query domain
	# Adjust the last residue based on the last residue index in the alignment and the last of the query domain
	#my $queryAdjustedStartIndex = MAX( $queryStartIndex, $startIndicesRef->[$QUERY_INDEX]);
	#my $queryAdjustedLastIndex  = MIN( $queryLastIndex,  $lastIndicesRef->[$QUERY_INDEX]); # get either the last residue index from the domain or the largest residue in the alignment
	my $queryAdjustedStartIndex = $queryStartIndex;
	if( $startIndicesRef->[$QUERY_INDEX] > $queryStartIndex){
	    $queryAdjustedStartIndex = $startIndicesRef->[$QUERY_INDEX];
	    DEBUG( "\t\tAdjusted query start index to $queryAdjustedStartIndex");
	    if( $queryAdjustedStartIndex < 0){
		die("\t\tqueryAdjustedStartIndex: $queryAdjustedStartIndex < 0");
	    }
	}
	# get either the last residue index from the domain or the largest residue in the alignment
	my $queryAdjustedLastIndex = $queryLastIndex;
	if( $lastIndicesRef->[$QUERY_INDEX] < $queryLastIndex){
	    $queryAdjustedLastIndex = $lastIndicesRef->[$QUERY_INDEX];
	    DEBUG( "\t\tAdjusted query last index to $queryAdjustedLastIndex");
	    if( $queryAdjustedLastIndex - $startIndicesRef->[$QUERY_INDEX] + 1 > $alignmentTotalLength ){
		die("\t\tqueryAdjustedLastIndex - startIndicesRef->[QUERY_INDEX] + 1: $queryAdjustedLastIndex - $startIndicesRef->[$QUERY_INDEX] + 1 (".($queryAdjustedLastIndex - $startIndicesRef->[$QUERY_INDEX] + 1).") > hsp->length('total'): ".$alignmentTotalLength);
	    }
	}

	
	# check if the domain is part of the alignment
	if( $queryAdjustedLastIndex < $queryAdjustedStartIndex){
	    DEBUG("\t\tDomain not part of alignment ( $queryAdjustedLastIndex < $queryAdjustedStartIndex) (query domain indices: $queryStartIndex..$queryLastIndex, start index: $startIndicesRef->[$QUERY_INDEX], last residue index in alignment: $lastIndicesRef->[$QUERY_INDEX])");
	    next;
	}

	# get subject indices corresponding with the adjusted query indices
	my $subjectStartIndex = $queryIndex2SubjectIndex[ $queryAdjustedStartIndex - $startIndicesRef->[$QUERY_INDEX]] + $queryIndicesAlignsToAGap[$queryAdjustedStartIndex - $startIndicesRef->[$QUERY_INDEX]];
	my $subjectLastIndex  = $queryIndex2SubjectIndex[ $queryAdjustedLastIndex - $startIndicesRef->[$QUERY_INDEX]];
	
	DEBUG("\t\tsubjectStartIndex: $subjectStartIndex; subjectLastIndex: $subjectLastIndex");
	
	if( $subjectStartIndex < 0){
	    die( "\t\tsubjectStartIndex: $subjectStartIndex < 0");
	}elsif( $subjectLastIndex- $startIndicesRef->[$SUBJECT_INDEX] + 1 > $alignmentTotalLength){
	    die( "\t\tsubjectLastIndex - startIndicesRef->[SUBJECT_INDEX] + 1: $subjectLastIndex - $startIndicesRef->[$SUBJECT_INDEX] + 1 (".($subjectLastIndex - $startIndicesRef->[$SUBJECT_INDEX] + 1).") > hsp->length('total'): ".$alignmentTotalLength);
	}

	for( my $subjectDomainIndex = 0; $subjectDomainIndex < $numberOfSubjectDomains; $subjectDomainIndex++){
	    DEBUG("\t\tsubjectDomainIndex: $subjectDomainIndex (".$domainIndicesRef->{ $labels[$SUBJECT_INDEX] }[$subjectDomainIndex][$DOMAIN_NAME_INDEX].": $domainIndicesRef->{ $labels[$SUBJECT_INDEX] }[$subjectDomainIndex][$START_INDEX]..$domainIndicesRef->{ $labels[$SUBJECT_INDEX] }[$subjectDomainIndex][$LAST_INDEX])");
	    if( $domainIndicesRef->{ $labels[$SUBJECT_INDEX] }[$subjectDomainIndex][$DOMAIN_NAME_INDEX] eq $queryDomainName){
		# my $subjectAdjustedStartIndex = MAX( $subjectStartIndex,  $domainIndicesRef->{ $labels[$SUBJECT_INDEX] }[$subjectDomainIndex][$START_INDEX]);
		# my $subjectAdjustedLastIndex  = MIN( $subjectLastIndex,   $domainIndicesRef->{ $labels[$SUBJECT_INDEX] }[$subjectDomainIndex][$LAST_INDEX]);
		
		my $subjectAdjustedStartIndex = $subjectStartIndex;
		if( $domainIndicesRef->{ $labels[$SUBJECT_INDEX] }[$subjectDomainIndex][$START_INDEX] > $subjectStartIndex){
		    $subjectAdjustedStartIndex = $domainIndicesRef->{ $labels[$SUBJECT_INDEX] }[$subjectDomainIndex][$START_INDEX];
		    DEBUG( "\t\t\tAdjusted subject start index to $subjectAdjustedStartIndex (for domain index $subjectDomainIndex)");
		    if( $subjectAdjustedStartIndex < 0){
			die("\t\t\tsubjectAdjustedStartIndex: $subjectAdjustedStartIndex < 0");
		    }elsif( $subjectAdjustedStartIndex < $startIndicesRef->[$SUBJECT_INDEX]){


			
			# FUTURE: Is this case actually valid?


			
			die("\t\t\tsubjectAdjustedStartIndex: $subjectAdjustedStartIndex < startIndicesRef->[SUBJECT_INDEX]: $startIndicesRef->[$SUBJECT_INDEX]");
		    }
		}
		my $subjectAdjustedLastIndex = $subjectLastIndex;
		if( $domainIndicesRef->{ $labels[$SUBJECT_INDEX] }[$subjectDomainIndex][$LAST_INDEX] < $subjectLastIndex){
		    $subjectAdjustedLastIndex = $domainIndicesRef->{ $labels[$SUBJECT_INDEX] }[$subjectDomainIndex][$LAST_INDEX];
		    DEBUG( "\t\t\tAdjusted subject last index to $subjectAdjustedLastIndex");
		    if( $subjectAdjustedLastIndex - $startIndicesRef->[$SUBJECT_INDEX] + 1 > $alignmentTotalLength ){
			die("\t\t\tsubjectAdjustedLastIndex - startIndicesRef->[SUBJECT_INDEX] + 1: $subjectAdjustedLastIndex - $startIndicesRef->[$SUBJECT_INDEX] + 1 (".($subjectAdjustedLastIndex - $startIndicesRef->[$SUBJECT_INDEX] + 1).") > hsp->length('total'): ".$alignmentTotalLength);
		    }elsif( $subjectAdjustedLastIndex > $lastIndicesRef->[$SUBJECT_INDEX]){


			
			# FUTURE: Is this case actually valid?


			
			die("\t\t\tsubjectAdjustedLastIndex: $subjectAdjustedLastIndex > lastIndicesRef->[SUBJECT_INDEX]: $lastIndicesRef->[$SUBJECT_INDEX]");
		    }
		}
		
		my $subjectDomainOverlap = $subjectAdjustedLastIndex - $subjectAdjustedStartIndex + 1;
		if( $subjectDomainOverlap > 0){
		    
		    my $queryAdjustedAdjustedStartIndex = $subjectIndex2QueryIndex[ $subjectAdjustedStartIndex - $startIndicesRef->[$SUBJECT_INDEX]] + $subjectIndicesAlignsToAGap[ $subjectAdjustedStartIndex - $startIndicesRef->[$SUBJECT_INDEX]];
		    if( $queryAdjustedAdjustedStartIndex != $queryAdjustedStartIndex){
			DEBUG("\t\t\tAdjusted query start index (for subject indices) to $queryAdjustedAdjustedStartIndex");
		    }
		    my $queryAdjustedAdjustedLastIndex = $subjectIndex2QueryIndex[ $subjectAdjustedLastIndex - $startIndicesRef->[$SUBJECT_INDEX]];
		    if( $queryAdjustedAdjustedLastIndex != $queryAdjustedLastIndex){
			DEBUG("\t\t\tAdjusted query last index (for subject indices) to $queryAdjustedAdjustedLastIndex");
		    }
		    my $queryDomainOverlapResidues = $queryAdjustedAdjustedLastIndex - $queryAdjustedAdjustedStartIndex + 1;
		    DEBUG("\t\t\tFound overlap: query: $queryAdjustedAdjustedStartIndex..$queryAdjustedAdjustedLastIndex ($queryDomainOverlapResidues); subject: $subjectAdjustedStartIndex..$subjectAdjustedLastIndex ($subjectDomainOverlap), domain index $subjectDomainIndex ($domainIndicesRef->{ $labels[$SUBJECT_INDEX] }[$subjectDomainIndex][$START_INDEX]..$domainIndicesRef->{ $labels[$SUBJECT_INDEX] }[$subjectDomainIndex][$LAST_INDEX])");

		    $queryOverlappingResidues += $queryDomainOverlapResidues;
		    
		    if( defined( $domainOverlapsRef)){
			# store overlap count into existing 2D array
			if( $queryDomainOverlapResidues > $domainOverlapsRef->[$queryDomainIndex][$subjectDomainIndex]){
			    if( $domainOverlapsRef->[$queryDomainIndex][$subjectDomainIndex] != 0){
				DEBUG("\t\t\tReplaced existing overlap for [$queryDomainIndex][$subjectDomainIndex]: $domainOverlapsRef->[$queryDomainIndex][$subjectDomainIndex]");
			    }
			    $domainOverlapsRef->[$queryDomainIndex][$subjectDomainIndex] = $queryDomainOverlapResidues;
			}else{
			    DEBUG( "\t\t\tSmaller than existing overlap: $domainOverlapsRef->[$queryDomainIndex][$subjectDomainIndex]");
			}
		    }
		}
	    }
	}# END foreach subject domain
    } # END foreach query domain

    my $overlap = ($queryOverlappingResidues / $queryTotalResidues) * 100.0;
    DEBUG( "\tHSP overlap: $overlap\%");
    if( $overlap < 0 ||
	$overlap > 100){
	die("HSP overlap: $overlap!");
    }
    
    return wantarray() ? ($queryOverlappingResidues, $queryTotalResidues) : $overlap;
} # END  calculateHspOverlap()


# calculate the overlap:
# foreach domain, for both query and subject
# + find the block index for the start of the domain
# + count the number of residues from the domain that are aligned
sub calculatePsiSemiGlobalOverlap{
    my ($taxon, $hitLabel, $startIndices, $lastIndices, $domainIndicesRef) = @_; 
    DEBUG( "calculatePsiSemiGlobalOverlap()");
    
    # first index is type
    my $domainStartIndices = [[]];
    my $domainLastIndices = [[]];
    # my $domainAlignmentStartIndices = [[]];
    # my $domainAlignmentLastIndices = [[]];

    my $numQueryDomainResidues = getQuerySubjectDomainIntersection( $domainIndicesRef, uc($taxon), uc($hitLabel), $domainStartIndices, $domainLastIndices);
    
    # Foreach domain
    # + Find the block index, i, for the start of the domain for the query
    # + While the domain region is part of the alignment for block i:
    # + For block i, calculate the overlap:
    #   * Determine if the domain for the subject sequence overlaps with block i
    #   * If it does, overlap += MAX(0,
    #                                MIN(query_last_alignment_index (for i), subject_last_alignment_index (for i)) - MAX( query_first_alignment_index (for i), subject_first_alignment_index (for i)) + 1)
    # + overlap /= query_domain_size
    
    if( $#{$domainStartIndices->[$QUERY_INDEX]} != $#{$domainLastIndices->[$QUERY_INDEX]}){
	die( "$#{$domainStartIndices->[$QUERY_INDEX]} != $#{$domainLastIndices->[$QUERY_INDEX]}");
    }
    DEBUG( "\tdomainStartIndices->[QUERY_INDEX]: @{$domainStartIndices->[$QUERY_INDEX]}");
    DEBUG( "\tdomainLastIndices->[QUERY_INDEX]:  @{$domainLastIndices->[$QUERY_INDEX]}");
    DEBUG( "\tdomainStartIndices->[SUBJECT_INDEX]: @{$domainStartIndices->[$SUBJECT_INDEX]}");
    DEBUG( "\tdomainLastIndices->[SUBJECT_INDEX]:  @{$domainLastIndices->[$SUBJECT_INDEX]}");
    
    my $overlap = 0.0;
    for( my $domainI = 0; $domainI <= $#{$domainStartIndices->[$QUERY_INDEX]}; $domainI++){
	DEBUG( "\t\tdomainI: $domainI");
	# + find the block index for the start of the domain
	my $blockI = 0;
	if( $domainLastIndices->[$QUERY_INDEX][$domainI] < $startIndices->[$QUERY_INDEX][$blockI]){
	    # domain ends before alignment starts 
	    DEBUG( "\t\tDomain ends at $domainLastIndices->[$QUERY_INDEX][$domainI], which is before the start of the alignment at index, $startIndices->[$QUERY_INDEX][$blockI]");
	    next;
	}
	
	#DEBUG("\t\twhile( $blockI <= $#{$lastIndices->[$QUERY_INDEX]} && $domainStartIndices->[$QUERY_INDEX][$domainI] > $lastIndices->[$QUERY_INDEX][$blockI])");
	while( $blockI <= $#{$lastIndices->[$QUERY_INDEX]} &&
	       $domainStartIndices->[$QUERY_INDEX][$domainI] > $lastIndices->[$QUERY_INDEX][$blockI]){
	    $blockI++;
	    #DEBUG("\t\twhile( $blockI <= $#{$lastIndices->[$QUERY_INDEX]} && $domainStartIndices->[$QUERY_INDEX][$domainI] > $lastIndices->[$QUERY_INDEX][$blockI])");
	}
	
	if( $blockI > $#{$startIndices->[$QUERY_INDEX]}){
	    # domain starts after alignment ends
	    DEBUG("\t\tDomain starts at index $domainLastIndices->[$QUERY_INDEX][$domainI], which is after the last alignment index of $lastIndices->[$QUERY_INDEX][$#{$lastIndices->[$QUERY_INDEX]}]");
	    next;
	}
	DEBUG( "\t\tDomain start index $domainStartIndices->[$QUERY_INDEX][$domainI] is in block index $blockI (indices $startIndices->[$QUERY_INDEX][$blockI]..$lastIndices->[$QUERY_INDEX][$blockI])");
	
	DEBUG( "\t\tSubject domain indices: $domainStartIndices->[$SUBJECT_INDEX][$domainI]..$domainLastIndices->[$SUBJECT_INDEX][$domainI]");

	# + While the query domain region is part of the alignment for block i:
	#   * Verify for both query and subject that 1) the start of the domain is before or equal to the end of the alignment for block $blockI
	#                                            2) the end of the domain is before or equal to the start index of the alignment for block $blockI
	DEBUG("\t\twhile( $blockI <= $#{$startIndices->[$QUERY_INDEX]} && $domainStartIndices->[$QUERY_INDEX][$domainI] <= $lastIndices->[$QUERY_INDEX][$blockI] && $domainLastIndices->[$QUERY_INDEX][$domainI] >= $startIndices->[$QUERY_INDEX][$blockI])");
	while( $blockI <= $#{$startIndices->[$QUERY_INDEX]} &&
	       $domainStartIndices->[$QUERY_INDEX][$domainI]   <= $lastIndices->[$QUERY_INDEX][$blockI] &&
	       $domainLastIndices->[$QUERY_INDEX][$domainI]    >= $startIndices->[$QUERY_INDEX][$blockI]){
	    DEBUG("\t\tDomain index $domainI is at least part of the alignment in block $blockI");
	    
	    # overlap += MAX(0, MIN(query_last_alignment_index (for i), subject_last_alignment_index (for i)) - MAX( query_first_alignment_index (for i), subject_first_alignment_index (for i)) + 1)
	    
	    # $startOffset is the number of residues in the block alignment before either of the domains start
	    my $queryStartOffset   = MAX(0, $domainStartIndices->[$QUERY_INDEX][$domainI]   - $startIndices->[$QUERY_INDEX][$blockI]);
	    my $subjectStartOffset = MAX(0, $domainStartIndices->[$SUBJECT_INDEX][$domainI] - $startIndices->[$SUBJECT_INDEX][$blockI]);
	    my $startOffset = MAX( $queryStartOffset, $subjectStartOffset);

	    my $blockSize = $lastIndices->[$QUERY_INDEX][$blockI] - $startIndices->[$QUERY_INDEX][$blockI] + 1;
	    # assert
	    if( $startOffset > $blockSize){
		DEBUG( "\t\t\tThe start offset for the alignment for domain index $domainI, $startOffset, is past the end of block index $blockI: $startOffset > $blockSize (queryStartOffset: $queryStartOffset, subjectStartOffset: $subjectStartOffset) (probably the subject domain hasn't started yet)");
		$blockI++;
		next;
	    }
	    
	    # $lastOffset is the number of residues in the block alignment before either of the domains last
	    my $queryLastOffset   = MAX(0, $lastIndices->[$QUERY_INDEX][$blockI]   - $domainLastIndices->[$QUERY_INDEX][$domainI]);
	    my $subjectLastOffset = MAX(0, $lastIndices->[$SUBJECT_INDEX][$blockI] - $domainLastIndices->[$SUBJECT_INDEX][$domainI]);

	    if( $subjectLastOffset > $blockSize){
		DEBUG("\t\t\tSubject domain is done (subjectLastOffset: $subjectLastOffset > blockSize: $blockSize)");
		last;
	    }

	    my $lastOffset = MAX( $queryLastOffset, $subjectLastOffset);
	    
	    # assert
	    if( $lastOffset > $blockSize){
		DEBUG("\t\t\tThe last offset (from the end of the block) for the alignment for domain index $domainI, $lastOffset, is larger than the size of block index $blockI: $lastOffset > $blockSize (queryLastOffset: $queryLastOffset, subjectLastOffset: $subjectLastOffset) (Hmm, the queryLastOffset should NOT be greater than the blockSize; an error?)");
		$blockI++;
		next;
	    }

	    # the overlap contribution for this block is the size of the block minus the start and last offsets
	    DEBUG("\t\t\tBlock overlap: $overlap += MAX(0, $blockSize - $startOffset - $lastOffset)");
	    $overlap += MAX(0, $blockSize - $startOffset - $lastOffset);

	    $blockI++;
	    if( $DEBUG){
		if( $blockI <= $#{$startIndices->[$QUERY_INDEX]}){
		    DEBUG("\t\t\twhile( $blockI <= $#{$startIndices->[$QUERY_INDEX]} && $domainStartIndices->[$QUERY_INDEX][$domainI] <= $lastIndices->[$QUERY_INDEX][$blockI] && $domainLastIndices->[$QUERY_INDEX][$domainI] >= $startIndices->[$QUERY_INDEX][$blockI])");
		}
	    }
	}
    }
    $overlap = ($overlap / $numQueryDomainResidues) * 100.0;
    DEBUG("\toverlap: $overlap\%");
    if( $overlap < 0 ||
	$overlap > 100){
	die("Invalid (PSI-SemiGLOBAL) overlap: $overlap!");
    }
    return $overlap;
} # END  calculatePsiSemiGlobalOverlap()

sub parseBlastHits{
    parseSearchIOHits(  { -format => "blast"}, @_);


}
    
# sub parseHmmerHits{
#     parseSearchIOHits(  { -format => "hmmer", -version => 3}, @_);
# }
    

sub parseSearchIOHits{
    my ($argsForNew, $data, $eValuesRef, $origIndicesRef, $hitsOutputFileName, $ignoreSelfHit, $storeOverlap, $allowMultipleSeqs, $allowMultipleHSPs, $combineHSPs) = @_;
    
    use Bio::SearchIO;
    use Bio::SeqIO;
    
    if( ! defined( $storeOverlap)){
	$storeOverlap = 0;
    }

    DEBUG( "parseSearchIOHits( argsForNew, data, eValuesRef, origIndicesRef, $hitsOutputFileName, $ignoreSelfHit, $storeOverlap, $allowMultipleSeqs, $allowMultipleHSPs, $combineHSPs)");
    # use BioPerl to parse BLAST hits information; store it in 2D array: $data (with E-value stored in @eValues)

    my $hitsIndex = 0;

    my %labels;  # used to detect duplicates
    
    my $in = new Bio::SearchIO( -file => "$hitsOutputFileName", %$argsForNew );
    while( my $result = $in->next_result ){
    
    	my $taxon = uc($result->query_name());
    	if( $DEBUG){
    	    print STDERR "DEBUG: taxon: $taxon\n";
    	    print STDERR "DEBUG: \thits: ".$result->num_hits()."\n";
    	}
	my %queryDomainNames;
	if( $storeOverlap){
	    # make of a hash of all of the query domain names (so that later we can do a quick check with the subject domain names for a match)
	    my $numberOfQueryDomains = $#{$domainIndices{ $taxon}} + 1;
	    for( my $queryDomainIndex = 0; $queryDomainIndex < $numberOfQueryDomains; $queryDomainIndex++){
		$queryDomainNames{ $domainIndices{ $taxon}[$queryDomainIndex][$DOMAIN_NAME_INDEX] } = 1;
	    }
	}
    
    	while( my $hit = $result->next_hit ){
    	    my $hitLabel = uc($hit->name);
    	    $hitLabel =~ s/^LCL\|//;  # strip of leading "lcl|" if present
    	    if( $ignoreSelfHit &&
    		$taxon eq $hitLabel){
    		if( $DEBUG){ print STDERR "DEBUG: Ignoring self hit \"$hitLabel\"\n"; }
    		# bioperl seems to have problems with self hits (the simpleAlign object I think, used to calculate the overlap)
    		next;
    	    }
	    if( $allowMultipleSeqs == $FALSE ){
		if( exists( $labels{$hitLabel})){
		    DEBUG( "Duplicate label: $hitLabel (skipping)");
		    next;
		}else{
		    $labels{ $hitLabel} = 1;
		}		    
	    }

	    if( $combineHSPs){
		# use best E-value
		my $bestEValue = "";
		my @hsps;
		while( my $hsp = $hit->next_hsp ){
		    my $eValue = $hsp->evalue;
		    $eValue =~ s/,$//;  # strip off the tailing comma
		    if( $bestEValue eq "" ||
			$eValue < $bestEValue){
			$bestEValue = $eValue;
		    }
		    push( @hsps, $hsp);
		}
		if( $#hsps < 0){
		    die( "ERROR: No HSPs found for $hitLabel!");
		}
		
		DEBUG( "\t$hitLabel\tbest evalue: $bestEValue");
		$data->[$hitsIndex][$TAXON_LABEL] = $hitLabel;
		#$data->[$hitsIndex][$EVALUE]     = $eValue;
		push( @$eValuesRef, $bestEValue);
		if( $storeOverlap){
		    $data->[$hitsIndex][$OVERLAP] = calculateMultiHspOverlap( $taxon, $hitLabel, \@hsps, \%domainIndices, \%queryDomainNames);
		}
		$hitsIndex++;
	    }else{
		
		while( my $hsp = $hit->next_hsp ){
		    my $eValue = $hsp->evalue;
		    $eValue =~ s/,$//;  # strip off the tailing comma
		    DEBUG( "\t$hitLabel\thsp evalue: $eValue");
		    
		    $data->[$hitsIndex][$TAXON_LABEL] = $hitLabel;
		    #$data->[$hitsIndex][$EVALUE]     = $evalue;
		    push( @$eValuesRef, $eValue);
		    if( $storeOverlap){
			$data->[$hitsIndex][$OVERLAP] = calculateHspOverlap( $taxon, $hitLabel, $hsp, \%domainIndices);
		    }
		    $hitsIndex++;
		    if( $allowMultipleHSPs == $FALSE){
			# skip if subsequent HSPs
			if( $DEBUG){
			    my @eValues;
			    while( my $hsp = $hit->next_hsp ){
				my $eValue = $hsp->evalue;
				push( @eValues, $eValue);
			    }
			    if( $#eValues >= 0){
				DEBUG( "\t\tSkipping over HSPs with the following E-values:" . join( ' ', @eValues));
			    }
			}
			last; # HSP
		    }
		} # END HSP loop
	    } # END if( combineHSPs)
    	} # END hit loop
    }
    
    DEBUG( "Done reading hits file");
    
    for( my $i = $#{$eValuesRef}; $i >= 0; $i--){
    	$origIndicesRef->[$i] = $i;
    }
    
    if( $DEBUG){ print STDERR "\@{\$eValuesRef} (\$#{\$eValuesRef}: $#{$eValuesRef}): \"@{$eValuesRef}\"\n"; }
    
    sort2Arrays($eValuesRef, $origIndicesRef);
    
    if( $DEBUG){ print STDERR "DEBUG: Done sorting ". ($#{$eValuesRef} + 1)." hits\n"; }
    if( $DEBUG){ print STDERR "DEBUG: Sorted original indices: @{$origIndicesRef}\n"; }
    if( $DEBUG){
    	print STDERR "\@{eValuesRef}: @{$eValuesRef}\n";
    	#print STDERR "Sorted: ";
    	#print "$eValues[$origIndices[0]] ";
    	for( my $i = 1; $i <= $#{$origIndicesRef}; $i++){
    	    #print STDERR "$eValues[$origIndices[$i]] ";
    	    if( $eValuesRef->[$i - 1] > $eValuesRef->[$i]){
    		die("ERROR: sorted array is NOT in order: Comparing ".($i-1)." and $i: $eValues[$i - 1] > $eValues[$i]; ");
    	    }
    	}
    	#print STDERR "\n";
    }
} # END parseSearchIOHits()


sub parsePsiSemiGlobalHits{
    my ($data, $eValuesRef, $origIndicesRef, $hitsOutputFileName, $ignoreSelfHit, $storeOverlap, $taxon, $allowMultipleSeqs) = @_;

    DEBUG( "parsePsiSemiGlobalHits( data, eValuesRef, origIndicesRef, $hitsOutputFileName, $ignoreSelfHit, $storeOverlap, $taxon, $allowMultipleSeqs)");
    
    # Parse the BLAST summary hits information; store it in 2D array: $data (with E-value stored in @eValues) and parse the 

    my $hitsIndex = 0;

    my $fileContents = fileToString( $hitsOutputFileName);

    # remove header (and any rounds that are not the last one)
    #$fileContents =~ s/^.+Sequences producing significant alignments[^\n]+\s*//m;
    $fileContents =~ s/^.+Sequences producing significant alignments[^\n]+\s*//s;  # "?" in case there's more than one round
    #if( $DEBUG){ print STDERR "fileContents: ".substr($fileContents, 0, 100)."\n"; }

    # copy the summary information from the contents
    my $summaryInfo = $fileContents;
    $summaryInfo =~ s/\n\n+.+$//s;
    
    my @fileContents = split( /\s*\n\s*/, $summaryInfo);
    if( $DEBUG){ print STDERR "Found ".($#fileContents + 1) . " summary lines in $hitsOutputFileName\n"; }

    # used to detect duplicates
    my %labels;
    my $NOT_SKIPPED = 0;
    my $SKIPPED = 1;
    my @skippedIndices;
    
    my $prevEValue = 0;
    for( my $lineI = 0; $lineI <= $#fileContents; $lineI++){
	# get taxon label and strip off leading "lcl|" if present
	# example: lcl|30749733  unnamed protein product                                  0.0    8e-62  8e-66    inf
        # example: lcl|d1u7pa_  c.108.1.17 (A:) Magnesium-dependent phosphatase-1,  ...   0.0    10545  1.00    9777	
	if( $fileContents[$lineI] =~ m/^(lcl\|)?(\S+)\s\s+(.+?)\s\s+(\d[\d\.]+)\s\s+(\S+)\s+(\S+)/){
	    my $hitLabel = uc($2);
	    my $eValue = $5;
	    
	    if( $DEBUG){ print STDERR "DEBUG:\t$hitLabel\thsp evalue: $eValue\n"; }
	    if( $allowMultipleSeqs == $FALSE ){
		if( exists( $labels{$hitLabel})){
		    DEBUG( "Duplicate label: $hitLabel (skipping)");
		    $skippedIndices[$hitsIndex] = $SKIPPED;
		    next;
		}else{
		    $skippedIndices[$hitsIndex] = $NOT_SKIPPED;
		    $labels{ $hitLabel} = 1;
		}		    
	    }
	    if( $eValue < $prevEValue){
		print STDERR "\nALERT:  E-value \"$eValue\" (from $fileContents[$lineI]) is lower than previous E-value \"$prevEValue\" (from ".$fileContents[$lineI - 1].")!\n\n";
	    }
	    $prevEValue = $eValue;
	    $data->[$hitsIndex][$TAXON_LABEL]           = $hitLabel;
	    push( @$eValuesRef, $eValue);
	    $hitsIndex++;
	}elsif($fileContents[$lineI] !~ m/^\s*$/){
	    die( "ERROR: Mal-formed line: $fileContents[$lineI] in $hitsOutputFileName (lineI: $lineI)");
	}# END hit loop
    }
    
    if( $storeOverlap){
	# partition the rest of the output by sequence
	my $alignmentInfo = $fileContents;
	$alignmentInfo =~ s/^.+?\n\n+//s;
	if( $DEBUG){ print STDERR "alignmentInfo: ".substr($alignmentInfo, 0, 100)."\n"; }
	my @blockAlignmentsSets = split( />/, $alignmentInfo);
	# The first index is what's before the split pattern (i.e., nothing)
	DEBUG("Found ".($#blockAlignmentsSets + 1 - 1)." blocks sets");
	for( my $blockSetI = 1; $blockSetI <= $#blockAlignmentsSets; $blockSetI++){
	    #if( $DEBUG){ print STDERR "blockAlignmentsSets[$blockSetI]: ".substr($blockAlignmentsSets[$blockSetI], 0, 100)."\n"; }
	    my $hitIndex = $blockSetI - 1;	    # ASSUME that there's not multiple HSPs per blocks alignment
	    if( $skippedIndices[ $hitIndex] == $SKIPPED){
		next;
	    }
	    my $hitLabel = $data->[$hitIndex][$TAXON_LABEL];
	    my @blockAlignmentsContents = split( /Block[^\n]+\n/, $blockAlignmentsSets[ $blockSetI]);
	    # The first index is what's before the split pattern
	    if( $blockAlignmentsContents[0] !~ m/\b$hitLabel\b/i){
		die( "ERROR: Expecting to find taxon label \"$hitLabel\" in header of alignment portion of the PSI-SemiGLOBAL output file but found: $blockAlignmentsContents[0]");
	    }
	    #print STDERR "Found ".($#blockAlignmentsContents + 1 - 1)." blocks\n";
	    
	    # store the aligned strings for each block (and for both query and subject) into a separate array cell
	    my @type2 = ('Query', 'Sbjct');
	    my $alignmentStrs = [[]];  # first index is the type (0: query, 1: subject); second index is block number
	    my $startIndices = [[]];  # first index is the type (0: query, 1: subject); second index is block number
	    my $lastIndices = [[]];   # first index is the type (0: query, 1: subject); second index is block number
	    #print STDERR "Assuming that the 1st index of blockAlignmentsContents is header info: $blockAlignmentsContents[0]\n";
	    for( my $blockI = 1; $blockI <= $#blockAlignmentsContents; $blockI++){
		my $pattern = '\s+(\d+)\s+[a-zA-Z\-]+\s+(\d+)';
		for( my $typeI = 0; $typeI <= $#TYPE; $typeI++){
		    if( $blockAlignmentsContents[$blockI] !~ m/$type2[$typeI]$pattern/){
			die( "Mal-formed block alignment ($type2[$typeI])");
		    }
		    $startIndices->[$typeI][$blockI - 1] = $1 - 1;
		    #$alignmentStrs->[$typeI][$blockI - 1] = $2;
		    $lastIndices->[$typeI][$blockI - 1] = $2 - 1;
		    ##DEBUG("Parsed startIndices->[$typeI][$blockI - 1]: $startIndices->[$typeI][$blockI - 1]; alignmentStrs->[$typeI][$blockI - 1]: $alignmentStrs->[$typeI][$blockI - 1]; lastIndices->[$typeI][$blockI - 1]: $lastIndices->[$typeI][$blockI - 1]");
		    #DEBUG("Parsed startIndices->[$typeI][$blockI - 1]: $startIndices->[$typeI][$blockI - 1]; lastIndices->[$typeI][$blockI - 1]: $lastIndices->[$typeI][$blockI - 1]");
		}
	    }
	    $data->[$hitIndex][$OVERLAP] = calculatePsiSemiGlobalOverlap( $taxon, $hitLabel, $startIndices, $lastIndices, \%domainIndices);
	}
    }


    
    if( $DEBUG){ print STDERR "DEBUG: Done reading hits file\n"; }


    for( my $i = $#{$eValuesRef}; $i >= 0; $i--){
	$origIndicesRef->[$i] = $i;
    }

    if( $DEBUG){ print STDERR "\@{\$eValuesRef} (\$#{\$eValuesRef}: $#{$eValuesRef}): \"@{$eValuesRef}\"\n"; }

    sort2Arrays($eValuesRef, $origIndicesRef);

    if( $DEBUG){ print STDERR "DEBUG: Done sorting ". ($#{$eValuesRef} + 1)." hits\n"; }
    if( $FALSE && $DEBUG){
	print STDERR "\@{eValuesRef}: @{$eValuesRef}\n";
	#print STDERR "Sorted: ";
	#print "$eValues[$origIndices[0]] ";
	for( my $i = 1; $i <= $#{$origIndicesRef}; $i++){
	    #print STDERR "$eValues[$origIndices[$i]] ";
	    if( $eValuesRef->[$i - 1] > $eValuesRef->[$i]){
		die("ERROR: sorted array is NOT in order: Comparing ".($i-1)." and $i: $eValues[$i - 1] > $eValues[$i]; ");
	    }
	}
	#print STDERR "\n";
    }
} # END parsePsiSemiGlobalHits()




sub parseSimpleHits{
    my ($data, $eValuesRef, $origIndicesRef, $hitsOutputFileName, $allowMultipleSeqs) = @_;

    DEBUG("parseSimpleHits($data, $eValuesRef, $origIndicesRef, $hitsOutputFileName, $allowMultipleSeqs)");
    
    # Parse just the summary BLAST hits information; store it in 2D array: $data (with E-value stored in @eValues)

    my $hitsIndex = 0;

    my $fileContents = fileToString( $hitsOutputFileName);

    # remove header
    #$fileContents =~ s/^.+Sequences producing significant alignments[^\n]+\s*//m;
    $fileContents =~ s/^.+?Sequences producing significant alignments[^\n]+\s*//s;  # "?" in case there's more than one round
    # remove detailed hits information
    # ASSUMES that the summary information has all of the hits (i.e., -num_descriptions is sufficently large)
    $fileContents =~ s/\n\n+.+$//s;
    
    my @fileContents = split( /\s*\n\s*/, $fileContents);
    if( $DEBUG){ print STDERR "Found ".($#fileContents + 1) . " summary lines in $hitsOutputFileName\n"; }

    # used to detect duplicates
    my %labels;
    
    my $prevEValue = 0;
    for( my $lineI = 0; $lineI <= $#fileContents; $lineI++){
	# get taxon label and strip off leading "lcl|" if present
	# example: lcl|30749733  unnamed protein product                                  0.0    8e-62  8e-66    inf
        # example: lcl|d1u7pa_  c.108.1.17 (A:) Magnesium-dependent phosphatase-1,  ...   0.0    10545  1.00    9777	
	if( $fileContents[$lineI] =~ m/^(lcl\|)?(\S+)\s\s+(.+?)\s\s+(\d[\d\.]+)\s\s+(\S+)/){
	    my $hitLabel = uc($2);
	    my $eValue = $5;
	    
	    if( $DEBUG){ print STDERR "DEBUG:\t$hitLabel\thsp evalue: $eValue\n"; }
	    if( $allowMultipleSeqs == $FALSE ){
		if( exists( $labels{$hitLabel})){
		    DEBUG( "Duplicate label: $hitLabel (skipping)");
		    next;
		}else{
		    $labels{ $hitLabel} = 1;
		}		    
	    }
	    if( $eValue < $prevEValue){
		print STDERR "\nALERT:  E-value \"$eValue\" (from $fileContents[$lineI]) is lower than previous E-value \"$prevEValue\" (from ".$fileContents[$lineI - 1].")!\n\n";
	    }
	    $prevEValue = $eValue;
	    $data->[$hitsIndex][$TAXON_LABEL]           = $hitLabel;
	    push( @$eValuesRef, $eValue);
	    $hitsIndex++;
	}elsif($fileContents[$lineI] !~ m/^\s*$/){
	    die( "ERROR: Mal-formed line: $fileContents[$lineI] in $hitsOutputFileName (lineI: $lineI)");
	    
	}# END hit loop
    }

    if( $DEBUG){ print STDERR "DEBUG: Done reading hits file\n"; }


    for( my $i = $#{$eValuesRef}; $i >= 0; $i--){
	$origIndicesRef->[$i] = $i;
    }

    if( $DEBUG){ print STDERR "\@{\$eValuesRef} (\$#{\$eValuesRef}: $#{$eValuesRef}): \"@{$eValuesRef}\"\n"; }

    sort2Arrays($eValuesRef, $origIndicesRef);

    if( $DEBUG){ print STDERR "DEBUG: Done sorting ". ($#{$eValuesRef} + 1)." hits\n"; }
    if( $FALSE && $DEBUG){
	print STDERR "\@{eValuesRef}: @{$eValuesRef}\n";
	#print STDERR "Sorted: ";
	#print "$eValues[$origIndices[0]] ";
	for( my $i = 1; $i <= $#{$origIndicesRef}; $i++){
	    #print STDERR "$eValues[$origIndices[$i]] ";
	    if( $eValuesRef->[$i - 1] > $eValuesRef->[$i]){
		die("ERROR: sorted array is NOT in order: Comparing ".($i-1)." and $i: $eValues[$i - 1] > $eValues[$i]; ");
	    }
	}
	#print STDERR "\n";
    }
} # END parseSimpleHits()


sub parseGlobalHits{
    my ($data, $eValuesRef, $origIndicesRef, $hitsOutputFileName, $totalSeqs, $allowMultipleSeqs) = @_;

    DEBUG( "parseGlobalHits($data, $eValuesRef, $origIndicesRef, $hitsOutputFileName, $totalSeqs, $allowMultipleSeqs)");
    
    # Parse just the hits (summary) information; store it in 2D array: $data (with E-value stored in @eValues)

    my $hitsIndex = 0;

    my $fileContents = fileToString( $hitsOutputFileName);

    # remove header
    $fileContents =~ s/^\d+\s+//s;
    my $headers = $fileContents;
    $headers =~ s/^([^\n]+)\n.+$/$1/s;
    $fileContents =~ s/^([^\n]+)\n//s;
    my @headers = split( /\t+/, $headers);
    if( $#headers < 0){
	die( "ERROR: Expecting third line of GLOBAL output file \"$hitsOutputFileName\" to be headers (found: @headers);");
    }
    my @headerStrs = (
	"query#",
	"gi-number",
	"score",
	"query length",
#	"P-value",
#	"E-value",
#	"hit_start_indices",
#	"probe_start_indices");
	"first AA of query",
	"last AA of query",
	"first AA of PSSM",
	"last AA of PSSM",
	"P-value",
	"P-value error");
    

    if( $#headerStrs != $#headers){
	print STDERR "ERROR:\nexpecting:\t".join(':',@headerStrs)."\nfound:\t\t".join(':',@headers)."\n";
	die( "ERROR: Differing numbers of header entries: expecting: $#headerStrs, found: $#headers;");
    }

    for( my $i = 0; $i <= $#headerStrs; $i++){
	if( $headerStrs[$i] ne $headers[$i]){
	    die( "ERROR: Header entries don't match: expecting: $headerStrs[$i], found: $headers[$i];");
	}
    }

    my @fileContents = split( /\s*\n\s*/, $fileContents);
    if( $DEBUG){ print STDERR "Found ".($#fileContents + 1) . " lines in $hitsOutputFileName\n"; }

    # used to detect duplicates
    my %labels;
    
    for( my $lineI = 0; $lineI <= $#fileContents; $lineI++){
	# 1	442592	122	129	1	129	1	151	0.99956	1e-50 in cd00001.hits
	if( $fileContents[$lineI] =~ m/^\d+\t(\S+)\t\d+\t\d+\t\d+\t\d+\t\d+\t\d+\t(\S+)\t(\S+)$/){
	    my $hitLabel = uc($1);
	    my $eValue = $2 * $totalSeqs;
	    if( $DEBUG){ print STDERR "DEBUG:\t$hitLabel\tE-value: $eValue\n"; }
	    if( $allowMultipleSeqs == $FALSE ){
		if( exists( $labels{$hitLabel})){
		    DEBUG( "Duplicate label: $hitLabel (skipping)");
		    next;
		}else{
		    $labels{ $hitLabel} = 1;
		}		    
	    }
	    $data->[$hitsIndex][$TAXON_LABEL]           = $hitLabel;
	    push( @$eValuesRef, $eValue);
	    $hitsIndex++;
	}elsif($fileContents[$lineI] !~ m/^\s*$/){
	    die( "ERROR: Mal-formed line: \"$fileContents[$lineI]\" in $hitsOutputFileName (lineI: $lineI)");
	    
	}# END hit loop
    }

    if( $DEBUG){ print STDERR "DEBUG: Done reading hits file\n"; }

    for( my $i = $#{$eValuesRef}; $i >= 0; $i--){
	$origIndicesRef->[$i] = $i;
    }

    if( $DEBUG){ print STDERR "\@{\$eValuesRef} (\$#{\$eValuesRef}: $#{$eValuesRef}): \"@{$eValuesRef}\"\n"; }
    
    sort2Arrays($eValuesRef, $origIndicesRef);

    if( $DEBUG){ print STDERR "DEBUG: Done sorting ". ($#{$eValuesRef} + 1)." hits\n"; }
    if( $FALSE && $DEBUG){
	print STDERR "\@{eValuesRef}: @{$eValuesRef}\n";
	#print STDERR "Sorted: ";
	#print "$eValues[$origIndices[0]] ";
	for( my $i = 1; $i <= $#{$origIndicesRef}; $i++){
	    #print STDERR "$eValues[$origIndices[$i]] ";
	    if( $eValuesRef->[$i - 1] > $eValuesRef->[$i]){
		die("ERROR: sorted array is NOT in order: Comparing ".($i-1)." and $i: $eValues[$i - 1] > $eValues[$i]; ");
	    }
	}
	#print STDERR "\n";
    }
} # END parseGlobalHits()


sub parseDomainLocs{
    # Put domain information into a 2D array (1 row per domains; columns are: domain name, start, and end)
    # Store the domain 2D array into a hash (keyed on the sequence label)
    # Convert to locations to indices
    my ($domainLocsFilename, $domainIndicesRef, $taxon) = @_;

    DEBUG( "parseDomainLocs( $domainLocsFilename, domainIndicesRef)");

    my @domainLocsLines = split( />/, fileToString($domainLocsFilename));
    for( my $lineI = 1; $lineI <= $#domainLocsLines; $lineI++){ # ignore index 0 (it's just what's before the first ">", right?)
	my @line = split(/\s+/, $domainLocsLines[$lineI]);
	if( $#line + 1 < 4){
	    die( "Mal-formed domain location line: @line (Doesn't even have a ".'taxon label\s+domainID\s+\d+\s+\d+)!'." (lineI: $lineI)\n");
	}
	my $taxonLabel = uc(shift @line);
	#print STDERR "DEBUG: taxonLabel: $taxonLabel\n";
	my $domainsArrayRef = [[]]; # 1 row per domains; columns are: domain name, start, and end
	if( ($#line + 1) % 3 != 0){
	    die( "Mal-formed domain location line: @line (Doesn't even fit the pattern:".'(domainID\s+\d+\s+\d+)+'." ($#line + 1: ".($#line + 1)."))!\n");
	}
	for( my $domainI = 0; $domainI <= $#line; $domainI += 3){
	    $domainsArrayRef->[$domainI / 3][$DOMAIN_NAME_INDEX] = uc($line[$domainI]);
	    $domainsArrayRef->[$domainI / 3][$START_INDEX]       = $line[$domainI + 1] - 1;  # "- 1" to convert to index
	    $domainsArrayRef->[$domainI / 3][$LAST_INDEX]        = $line[$domainI + 2] - 1;  # "- 1" to convert to index
	}
	$domainIndicesRef->{$taxonLabel} = $domainsArrayRef
    }
    DEBUG( "Found ".scalar keys (%{$domainIndicesRef}) ." entries in $domainLocsFilename");

}  # END parseDomainLocs()

# #
# # for tabular HMMer output
# #
# sub parseHmmerHits{
#     my ($data, $eValuesRef, $origIndicesRef, $hitsOutputFileName, $totalSeqs, $calculateOverlaps, $allowMultipleSeqs) = @_;
# 
#     DEBUG( "parseHmmerHits($data, $eValuesRef, $origIndicesRef, $hitsOutputFileName, $totalSeqs, $calculateOverlaps, $allowMultipleSeqs)");
#     
#     # Parse just the hits information; store it in 2D array: $data (with E-value stored in @eValues)
# 
#     my $hitsIndex = 0;
# 
#     my $fileContents = fileToString( $hitsOutputFileName);
# 
#     if( $fileContents !~ m/#\s*(target[^\n]+)/){
# 	die( "ERROR: Could not find header row in HMMER table output file \"$hitsOutputFileName\";");
#     }
#     my $headers = $1;    
#     my @headers = split( /\s\s+/, $headers);
#     if( $#headers < 0){
# 	die( "ERROR: Expecting header row in HMMER table output file \"$hitsOutputFileName\" (found: @headers);");
#     }
#     my @headerStrs = (
# 	"target name",
#         "accession",
# 	"query name",
# 	"accession",
# 	"E-value",
# 	"score",
# 	"bias"
# 	# E-value  score  bias   exp reg clu  ov env dom rep inc description of target
# 	);
# 
#     # if( $#headerStrs != $#headers){
#     # 	print STDERR "ERROR:\nexpecting:\t".join(':',@headerStrs)."\nfound:\t\t".join(':',@headers)."\n";
#     # 	die( "ERROR: Differing numbers of header entries: expecting: $#headerStrs, found: $#headers;");
#     # }
# 
#     if( $#headers < $#headerStrs){
# 	print STDERR "ERROR:\nexpecting (at least):\t".join(':',@headerStrs)."\nfound:\t\t".join(':',@headers)."\n";
# 	die( "ERROR: Differing numbers of header entries: expecting: $#headerStrs, found: $#headers;");
#     }
# 
#     for( my $i = 0; $i <= $#headerStrs; $i++){
# 	if( $headerStrs[$i] ne $headers[$i]){
# 	    die( "ERROR: Header entries don't match: expecting: $headerStrs[$i], found: $headers[$i];");
# 	}
#     }
# 	
#     # remove header lines
#     $fileContents =~ s/#.*\n//g;
#     
#     my @fileContents = split( /\s*\n\s*/, $fileContents);
#     if( $DEBUG){ print STDERR "Found ".($#fileContents + 1) . " lines in $hitsOutputFileName\n"; }
# 
#     # used to detect duplicates
#     my %labels;
#     
#     for( my $lineI = 0; $lineI <= $#fileContents; $lineI++){
# 	# 30749733             -          cd00001              -            1.8e-58  196.2   0.3   2.1e-58  196.1   0.2   1.0   1   0   0   1   1   1   1 -
# 	if( $fileContents[$lineI] =~ m/^(\S+)\s+\S+\s+\S+\s+\S+\s+(\S+)\s+/){
# 	    my $hitLabel = uc($1);
# 	    my $eValue = $2;
# 	    if( $DEBUG){ print STDERR "DEBUG:\t$hitLabel\tE-value: $eValue\n"; }
# 	    if( $allowMultipleSeqs == $FALSE ){
# 		if( exists( $labels{$hitLabel})){
# 		    DEBUG( "Duplicate label: $hitLabel (skipping)");
# 		    next;
# 		}else{
# 		    $labels{ $hitLabel} = 1;
# 		}		    
# 	    }
# 	    $data->[$hitsIndex][$TAXON_LABEL]           = $hitLabel;
# 	    push( @$eValuesRef, $eValue);
# 	    $hitsIndex++;
# 	}elsif($fileContents[$lineI] !~ m/^\s*$/){
# 	    die( "ERROR: Mal-formed line: \"$fileContents[$lineI]\" in $hitsOutputFileName (lineI: $lineI)");
# 	    
# 	}# END hit loop
#     }
# 
#     if( $DEBUG){ print STDERR "DEBUG: Done reading hits file\n"; }
# 
#     for( my $i = $#{$eValuesRef}; $i >= 0; $i--){
# 	$origIndicesRef->[$i] = $i;
#     }
# 
#     if( $DEBUG){ print STDERR "\@{\$eValuesRef} (\$#{\$eValuesRef}: $#{$eValuesRef}): \"@{$eValuesRef}\"\n"; }
# 
#     sort2Arrays($eValuesRef, $origIndicesRef);
# 
#     if( $DEBUG){ print STDERR "DEBUG: Done sorting ". ($#{$eValuesRef} + 1)." hits\n"; }
#     if( $FALSE && $DEBUG){
# 	print STDERR "\@{eValuesRef}: @{$eValuesRef}\n";
# 	#print STDERR "Sorted: ";
# 	#print "$eValues[$origIndices[0]] ";
# 	for( my $i = 1; $i <= $#{$origIndicesRef}; $i++){
# 	    #print STDERR "$eValues[$origIndices[$i]] ";
# 	    if( $eValuesRef->[$i - 1] > $eValuesRef->[$i]){
# 		die("ERROR: sorted array is NOT in order: Comparing ".($i-1)." and $i: $eValues[$i - 1] > $eValues[$i]; ");
# 	    }
# 	}
# 	#print STDERR "\n";
#     }
# } # END parseHmmerHits()


#
# Parser for HMMer output (including alignments)
#
# Populate data structure to mirror BioPerl SearchIO data structures by parsing hitsOutputFileName (except, combine all domains into a single alignment).
# eValuesRef and origIndicesRef are populated by parseSearchIOHits().
# If storeOverlap == 1, then populate the $OVERLAP column for each hit in data (again, via parseSearchIOHits()).

# Options and their affects:
# combineHSPs:        Ignored
# ignoreMultipleHSPs: Ignored
# allowMultipleHSPs:  Ignored
# ignoreMultipleSeqs: Honored
# allowMultipleSeqs:  Honored

sub parseHmmerHits{
    my ($data, $eValuesRef, $origIndicesRef, $hitsOutputFileName, $ignoreSelfHit, $storeOverlap, $allowMultipleSeqs, $allowMultipleHSPs, $combineHSPs, $taxon) = @_;

    DEBUG( "parseHmmerHits(data, eValuesRef, origIndicesRef, $hitsOutputFileName, $ignoreSelfHit, $storeOverlap, $allowMultipleSeqs, $allowMultipleHSPs, $combineHSPs, $taxon)");
    
    # Outline:
    # Split output on "\n>>"
    # Look at index 0 from split
    # Skip headers
    # # if ! allowMultipleHSPs (so, only consider the best domain/HSP)
    #   Record "full sequence" E-value
    #   Record label (to verify that it matches with the detailed output later)
    #   # Record number of domains
    #
    # If storeOverlap
    #   foreach hit (indices 1.. from split)
    #     Split on "\s+== domain\s*"
    #     Verify name (with summary info above)
    #     Record summary information
    #       Skip headers
    #       For each domain:
    #         # if ! allowMultipleHSPs
    #         #   Skip if not most significant domain/HSP
    #         # Record i-Evalue ("independent E-value: the significance of the sequence in the whole database search, if this were the only domain we had identified.")
    #         
    #         # Record from and to for the query sequence (hmm) and the subject sequence (ali)
    #      If storeOverlap, then process alignment information
    #        Combine all query and subject alignments together respectively
    #        Get total length
    #        Calculate overlap

    my $fileContents = fileToString( $hitsOutputFileName);
    my @hitsContents = split( /\n>>/, $fileContents);
    DEBUG( "Found $#hitsContents detailed (i.e., with alignments) hits in $hitsOutputFileName");
    
    #
    # Pull out E-values and labels from the summary information
    #
    # used to detect duplicates
    my %labels;
    my $NOT_SKIPPED = 0;
    my $SKIPPED = 1;
    my @skippedIndices;
    
    my $summaryContents = $hitsContents[0];
    $summaryContents =~ s/^.+\sDescription[\s\-]+//s;  # remove header info (from the beginning)
    $summaryContents =~ s/\s*\-+\s*inclusion\s+threshold.+$//s;  # remove insignificant hits (from the end)
    $summaryContents =~ s/\s*Domain annotation for each sequence \(and alignments\)\:\s*$//s;  # remove header to next section (from the end)
    DEBUG( "summaryContents has " . length( $summaryContents) . " chars: " . substr( $summaryContents, 0, 33) . "..." . substr($summaryContents, -33));

    my @hits = split( /\n/, $summaryContents);
    DEBUG( "Found $#hits significant (summary) hits");

    # --- full sequence ---   --- best 1 domain ---    -#dom-
    #  E-value  score  bias    E-value  score  bias    exp  N  Sequence                   Description
    #  ------- ------ -----    ------- ------ -----   ---- --  --------                   -----------
    # 1.6e-254  846.7  12.4   1.8e-254  846.5  12.4    1.0  1  up_Q8K367_Q8K367_MOUSE      
    # 6.8e-252  838.0  13.8   8.5e-252  837.7  13.8    1.1  1  pfam21_Q3TUT6_Q3TUT6_MOUSE  
    my $prevEValue = 0;
    my $hitsIndex = 0;
    foreach my $line (@hits){
	if( $line !~ m/^\s*([\deE\-\+\.]+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(.+)/){
	    die( "ERROR: Mal-formed summary line in $hitsOutputFileName: $line");
	}
	my $eValue = $1;
	my $hitLabel = uc($2);
	$hitLabel =~ s/\s+$//;  # remove trailing whitespace

	DEBUG( "\t$hitLabel\tE-value: $eValue");
	if( $allowMultipleSeqs == $FALSE ){
	    if( exists( $labels{$hitLabel})){
		die( "Duplicate label: $hitLabel!");
		DEBUG( "Duplicate label: $hitLabel (skipping)");
		$skippedIndices[$hitsIndex] = $SKIPPED;
		next;
	    }else{
		$skippedIndices[$hitsIndex] = $NOT_SKIPPED;
		$labels{ $hitLabel} = 1;
	    }		    
	}
	if( $eValue < $prevEValue){
	    print STDERR "\nALERT:  E-value \"$eValue\" (from $line) is lower than previous E-value \"$prevEValue\"!\n\n";
	}
	$prevEValue = $eValue;
	$data->[$hitsIndex][$TAXON_LABEL] = $hitLabel;
	push( @$eValuesRef, $eValue);
	$hitsIndex++;
    }
    DEBUG( "Done reading hits file");

    for( my $i = $#{$eValuesRef}; $i >= 0; $i--){
    	$origIndicesRef->[$i] = $i;
    }

    if( $DEBUG){ print STDERR "\@{\$eValuesRef} (\$#{\$eValuesRef}: $#{$eValuesRef}): \"@{$eValuesRef}\"\n"; }
    
    sort2Arrays($eValuesRef, $origIndicesRef);
    
    if( $DEBUG){ print STDERR "DEBUG: Done sorting ". ($#{$eValuesRef} + 1)." hits\n"; }
    if( $DEBUG){ print STDERR "DEBUG: Sorted original indices: @{$origIndicesRef}\n"; }
    if( $DEBUG){
    	print STDERR "\@{eValuesRef}: @{$eValuesRef}\n";
    	#print STDERR "Sorted: ";
    	#print "$eValues[$origIndices[0]] ";
    	for( my $i = 1; $i <= $#{$origIndicesRef}; $i++){
    	    #print STDERR "$eValues[$origIndices[$i]] ";
    	    if( $eValuesRef->[$i - 1] > $eValuesRef->[$i]){
    		die("ERROR: sorted array is NOT in order: Comparing ".($i-1)." and $i: $eValues[$i - 1] > $eValues[$i]; ");
    	    }
    	}
    	#print STDERR "\n";
    }
    
    if( $storeOverlap){
	#
	# Process each hit
	#
	my $numSignificantHits = $hitsIndex; # no "- 1" because we're converting from index to a count
	if( $#hitsContents + 1 < $numSignificantHits){
	    die("ERROR: Found $numSignificantHits in the summary info but only found " . ($#hitsContents + 1) . " detailed (i.e., alignment) hits!");
	}
	$hitsContents[ $#hitsContents] =~ s/\s*Internal pipeline.+$//s; # remove footer (in the event that there's no non-significant hits)
	    
	for( my $hitsContentsI = 1; $hitsContentsI <= $numSignificantHits; $hitsContentsI++){  # start with 1 to start to ignore the header & summary info
	    my $hitIndex = $hitsContentsI - 1;
	    if( $skippedIndices[ $hitIndex] == $SKIPPED){
		die( "Not sure what to do here!");
		next;
	    }
	    my $hitLabel = $data->[$hitIndex][$TAXON_LABEL];
	    
	    my @domains = split( /\s+==\s*domain\s*/, $hitsContents[ $hitsContentsI]);
	    my $domainI = 0;
	    
	    # The first index is what's before the split pattern
	    my $domainHeaderContents = $domains[ $domainI];
	    
	    if( $domainHeaderContents !~ m/^\s*$hitLabel\b/i){
		die( "ERROR: Expecting to find taxon label \"$hitLabel\" in header of alignment portion of the HMMer output file but found: $domainHeaderContents");
	    }
	    #
	    # get from and tos
	    #
	    # >> up_Q4RYL0_Q4RYL0_TETNG  
	    #    #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
	    #  ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
	    #    1 !  266.7   1.1   1.5e-80   3.4e-78       4     250 ..      62     308 ..      59     355 .. 0.93
	    #    2 !   15.0   0.1   0.00048      0.11     325     354 ..     429     458 .]     401     458 .] 0.82
	    # 
	    #   Alignments for each domain:
	    
	    $domainHeaderContents =~ s/^.+\n[\s\-]+\n//s;  # remove header information from the beginning
	    $domainHeaderContents =~ s/\s*Alignments for each domain:\s*$//s;  # remove footer to the header

	    # Arrays of array references. The first index is start or end index, second index is the TYPE {query, subject}
	    my @arrayOfStartIndices;
	    my @arrayOfLastIndices;
	    
	    my @domainHeaderSummaries = split( /\n/, $domainHeaderContents);
	    if( $#domainHeaderSummaries == 1 &&
		$domainHeaderSummaries[1] =~ m/No individual domains that satisfy reporting thresholds/i){
		DEBUG("Found a hit with no significant domains: $hitLabel");
		$#domainHeaderSummaries = -1;
		$#domains = -1;
	    }
	    # Future: record the number of domains from the initial summary information section and verify that it matches with the number found here
	    for( my $domainHeaderSummariesI = 0; $domainHeaderSummariesI <= $#domainHeaderSummaries; $domainHeaderSummariesI++){
		my $domainHeaderSummaryLine = $domainHeaderSummaries[ $domainHeaderSummariesI];
		if( $domainHeaderSummaryLine !~ m/^\s*\d+\s+([\!\?])\s+\S+\s+\S+\s+\S+\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\d+)\s+(\d+)/){
		    die("ERROR: Mal-formed domain summary line: $domainHeaderSummaryLine (domainHeaderSummariesI: $domainHeaderSummariesI (out of $#domainHeaderSummaries))");
		}
		if( $1 ne "!"){
		    DEBUG( "\tDomain index $domainHeaderSummariesI for $hitLabel is not significant (symbol: $1) (skipping)");
		    next;
		}
		my $startIndices = [];
		my $lastIndices = [];
		$startIndices->[$QUERY_INDEX]   = $2 - 1;
		$lastIndices->[$QUERY_INDEX]    = $3 - 1;
		$startIndices->[$SUBJECT_INDEX] = $4 - 1;
		$lastIndices->[$SUBJECT_INDEX]  = $5 - 1;
		
		$arrayOfStartIndices[ $domainHeaderSummariesI] = $startIndices;
		$arrayOfLastIndices[  $domainHeaderSummariesI] = $lastIndices;
		
		# for( my $typeI = 0; $typeI <= $#TYPE; $typeI++){
		#     DEBUG( "\tDomain index $domainHeaderSummariesI: startIndices[$TYPE[$typeI]]: $startIndices->[$typeI][$domainHeaderSummariesI], lastIndices[$TYPE[$typeI]]: $lastIndices->[$typeI][$domainHeaderSummariesI]");
		# }
	    }

	    my $aggregateOverlapCount = 0;
	    my $aggregateResidueCount = 0;	

	    for( $domainI = 1; $domainI <= $#domains; $domainI++){
		if( !defined( $arrayOfStartIndices[$domainI - 1]->[$QUERY_INDEX])){
		    #DEBUG("\tDomain $domainI is not significant");
		    next;
		}
		my @domainContents = split( /\n/, $domains[ $domainI]);
		my $lineI = 1;  # skip over first line
		# Expecting: (NOTE: HMMer may output a line before the query line (NYI))
		# 1 query line
		# 2 midline (ignore)
		# 3 subject consensus line
		# 4 ignore
		# 5 whitespace
		my $queryAlignment = '';
		my $subjectAlignment = '';
		while( $lineI + 3 <= $#domainContents){
		    my $CS_RFPattern = '^\s+\S+\s+(CS|RF)\s*$';
		    if( $lineI <= $#domainContents && $domainContents[$lineI] =~ m/$CS_RFPattern/){
			# CS / RF annotation lines are ignored
			$lineI++;
			if( $lineI + 3 > $#domainContents){
			    die( "ERROR: After adjusting for CS or RF annotation line, there's not enough lines to have a valid alignment of both the query and the subject sequences! (lineI: $lineI; largest index: $#domainContents)");
			}
		    }
		    #                                vvvvvvv iteration tag for the final iteration output
		    my $pattern = '^\s*' . $taxon . '(\-i\d+)?\s+(\d+|\-)\s+([a-zA-Z\.\-]+)\s+(\d+|\-)';
		    if( $domainContents[$lineI] !~ m/$pattern/i){
			die("ERROR: Mal-formed query alignment line: \"$domainContents[$lineI]\" did not match the pattern: $pattern (lineI: $lineI [out of $#domainContents]; already found queryAlignment: $queryAlignment)");
		    }
		    $queryAlignment .= $3;
		    $lineI++; # done wity query line
		    $lineI++; # skip over midline
		    $pattern = '^\s*' . $hitLabel . '\s+(\d+|\-)\s+([a-zA-Z\.\-]+)\s+(\d+|\-)';
		    if( $domainContents[$lineI] !~ m/$pattern/i){
			die("ERROR: Mal-formed subject alignment line: \"$domainContents[$lineI]\" did not match the pattern: $pattern (lineI: $lineI [out of $#domainContents]; already found subjectAlignment: $subjectAlignment)");
		    }
		    $subjectAlignment .= $2;
		    $lineI++; # done wity subject sequence line
		    $lineI++; # skip over HMMer posterior probabilities 
		    if( $lineI <= $#domainContents &&
			$domainContents[$lineI] !~ m/^\s*$/){
			die("ERROR: Expected empty line at line index $lineI in the alignment section for domain $hitsContentsI for $hitLabel, but found \"$domainContents[$lineI]\"");
		    }
		    $lineI++; # skip over blank line
		}
		# HMMer puts "."s into the query sequence to indicate a gap
		$queryAlignment =~ s/\./-/g;
		
		DEBUG("\tquery alignment: $queryAlignment ($taxon)");
		DEBUG("\tsbjct alignment: $subjectAlignment ($hitLabel)");


		my $alignmentTotalLength = length( $queryAlignment);
		my @alignmentSeqs = ($queryAlignment, $subjectAlignment);
		
		# calculate overlap for domain
		my ($overlapCount, $residueCount) = calculateAlignmentOverlap( $taxon, $hitLabel, $alignmentTotalLength, \@alignmentSeqs, $arrayOfStartIndices[ $domainI - 1], $arrayOfLastIndices[ $domainI - 1], \%domainIndices);
		$aggregateOverlapCount += $overlapCount;
		$aggregateResidueCount += $residueCount;
		DEBUG( "\taggregateOverlapCount: $aggregateOverlapCount; aggregateResidueCount: $aggregateResidueCount; (".($aggregateOverlapCount / $aggregateResidueCount * 100.0)."%)\n");
	    }
	    if( $aggregateResidueCount == 0){
		$data->[$hitIndex][$OVERLAP] = 0;
	    }else{
		$data->[$hitIndex][$OVERLAP] = $aggregateOverlapCount / $aggregateResidueCount * 100.0;
	    }
	} # END for( $hitsContentsI )
    } # END if( $storeOverlap) 
} # END parseHmmerHits()



sub fileToString{
    my $file_ = $_[0];  # Input file
    
    my $msg = "Input file '" . $file_ . "' failed to open.\n" . "Died";

    my $terminator = $/;
    undef $/;
    open(INPUT, "<$file_") or die($msg);
    my $str = <INPUT>; # terminator undefined : $str is the whole file.
    close INPUT;

    $/ = $terminator; # Restore for normal behavior later

    return $str;
}
    
