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
# Copyright 2016-2017 Hyrum D. Carroll

# Description: For each Uniprot accession (or Uniprot accession and species name tuple) output it's taxonomy

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

my $TEMP_DIR=".";  # Default to the current directory
my $NCBI_TAXONOMY_FTP_URL = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/";
my $NAMES_FILENAME_DEFAULT = "names.dmp";  # From NCBI Taxononmy
my $NODES_FILENAME_DEFAULT = "nodes.dmp";  # From NCBI Taxononmy
my $SPECIES_LIST_FILENAME_DEFAULT = "speclist.txt";  # From https://www.uniprot.org/docs/speclist.txt
my $SPECIES_LIST_COL_CODE = 0;   # column index for the species code
my $SPECIES_LIST_COL_TAXON = 2;  # column index for the taxon id
my $ACCESSION_2_TAXON_ID_FULL_URL = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping";  
my $ACCESSION_2_TAXON_ID_FULL_FILENAME = "idmapping_selected.tab";  # Original from $ACCESSION_2_TAXON_ID_FULL_URL/$ACCESSION_2_TAXON_ID_FULL_FILENAME # ~12 GB
my $ACCESSION_2_TAXON_ID_FILENAME_DEFAULT = "idmapping_selected-cached.tab";  # Parsed down output from uniprotId2taxonomy.pl (greatly speeds up repeated runs) (remove this file if new accession and/or species codes are required (in the uniprotTuplesFilename))

my $accession2taxonFilename; # default set later (with tempDir) if not set by the user
my $namesFilename;           # default set later (with tempDir) if not set by the user
my $nodesFilename;           # default set later (with tempDir) if not set by the user
my $speciesListFilename;     # default set later (with tempDir) if not set by the user
my $speciesListOldFilename = ""; # Older version of speclist (to get entries that have been deleted)
my $tempDir = $TEMP_DIR;
my $uniprotTuplesFilename = "";
my $uniprotSpeciesCode = "";
my $outputFilename = "";
my $outputSeparator = "\t";
my $man = 0;
my $help = 0;

unless( GetOptions("-v|verbose+"        => \$VERBOSE,
		   "-d|debug+"          => \$DEBUG,
		   "uniprots=s"         => \$uniprotTuplesFilename,
		   "speclist=s"         => \$speciesListFilename,
		   "speclist2=s"        => \$speciesListOldFilename,
		   "names=s"            => \$namesFilename,
		   "nodes=s"            => \$nodesFilename,
		   "uniprotMappings=s"  => \$accession2taxonFilename,
		   "tempDir=s"          => \$tempDir,
		   "out=s"              => \$outputFilename,
		   "separator=s"        => \$outputSeparator,
		   'help|?'             => \$help,
		   'man'                => \$man)){
    pod2usage(2);    
}

# set filename defaults (not done above because they depend on $tempDir)
if( ! defined( $accession2taxonFilename)){
    $accession2taxonFilename = "$tempDir/$ACCESSION_2_TAXON_ID_FILENAME_DEFAULT";
}

if( ! defined( $namesFilename) ){
    $namesFilename = "$tempDir/$NAMES_FILENAME_DEFAULT";
}

if( ! defined( $nodesFilename) ){
    $nodesFilename = "$tempDir/$NODES_FILENAME_DEFAULT";
}

if( ! defined( $speciesListFilename)){
    $speciesListFilename = "$tempDir/$SPECIES_LIST_FILENAME_DEFAULT";
}


DEBUG( "$0 parameters:
    uniprotTuplesFilename:        $uniprotTuplesFilename
    speciesListFilename:          $speciesListFilename
    speciesListOldFilename:       $speciesListOldFilename
    accession2taxonFilename:      $accession2taxonFilename
    namesFilename:                $namesFilename
    nodesFilename:                $nodesFilename
    tempDir:                      $tempDir
    outputFilename:               $outputFilename
    outputSeparator (\"s added):   \"$outputSeparator\"
    help:                         $help
    man:                          $man
    VERBOSE:                      $VERBOSE
    DEBUG:                        $DEBUG
");

pod2usage(-exitstatus => 0, -verbose => 1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

if( $uniprotTuplesFilename eq ""){
    pod2usage( -message => "ERROR: uniprotTuplesFilename not specified");
}
if( ! -s "$uniprotTuplesFilename" ){
    pod2usage( -message => "\nERROR: uniprotTuplesFilename: \"$uniprotTuplesFilename\" does not exist!\n");
}

if( ! -s "$speciesListFilename" ){
    # try to download the file from uniprot.org
    VERBOSE( "$speciesListFilename not found, downloading the species list from uniprot.org . . .");
    downloadFile( "https://www.uniprot.org/docs/speclist.txt", "$speciesListFilename");
}

if( ! -s "$speciesListFilename"){
    pod2usage( -message => "\nERROR: speciesListFilename: \"$speciesListFilename\" does not exist!\n");
}

my $usingSimplifiedAccession2taxon = $TRUE;  # Flag to indicate format of the accession 2 taxon mapping file.  When set, the format is <accession>\s<taxon ID> (see $accession2taxonRegex).  Otherwise, <accession>\t[11 other columns]\t<taxon ID>...  Also, if not set, produce a file of the subset used in the simplified format
my $accession2taxonRegex = qr/^(\S+)\t(\S+)$/;
if( ! -s "$accession2taxonFilename"){
    VERBOSE( "accession2taxonFilename \"$accession2taxonFilename\" does not exist.  Using $ACCESSION_2_TAXON_ID_FULL_FILENAME instead.");
    $accession2taxonFilename = "$tempDir/$ACCESSION_2_TAXON_ID_FULL_FILENAME";
    if( ! -s "$accession2taxonFilename"){
	VERBOSE( "$accession2taxonFilename not found."); #, downloading the accession to taxon mappings from $ACCESSION_2_TAXON_ID_FULL_URL . . .");
	print("Downloading the accession to taxon mappings from $ACCESSION_2_TAXON_ID_FULL_URL.  This may take a couple of minutes . . .");
	downloadFile( "$ACCESSION_2_TAXON_ID_FULL_URL/$ACCESSION_2_TAXON_ID_FULL_FILENAME.gz", "$accession2taxonFilename.gz");
	# if( -s "$accession2taxonFilename.gz"){ # handled by downloadFile()
	print("Decompressing the accession to taxon mappings file \"$accession2taxonFilename.gz\".  This may take a couple of minutes . . .");
	system( "gunzip $accession2taxonFilename.gz");
        #  }
    }
    $usingSimplifiedAccession2taxon = $FALSE;
    # get the 1st and 13th columns:
    $accession2taxonRegex = qr/^([^\t]+)\t(?:[^\t]*\t){11}([^\t]+)/; # "(?:regex)" is a non-capturing group (so only $1 and $2 will be set)
}

if( ! -s "$namesFilename" ||
    ! -s "$nodesFilename" ){
    my $taxTarBallBaseFilename = "taxdump.tar.gz";
    my $taxTarBallFilename = "$tempDir/$taxTarBallBaseFilename";
    if( ! -s "$taxTarBallFilename"){
	VERBOSE( "$namesFilename and/or $nodesFilename are not found, downloading them from $NCBI_TAXONOMY_FTP_URL . . .");
	downloadFile( "$NCBI_TAXONOMY_FTP_URL$taxTarBallBaseFilename", "$taxTarBallFilename");
    }
    # untar the tarball in the tempDir directory, then clean-up unused files
    my $cmd="tar -C $tempDir -xzf $taxTarBallFilename && rm $tempDir/{{citations,division,gencode,merged,delnodes}.dmp,gc.prt,readme.txt}";
    if( $DEBUG){ print "cmd: $cmd\n"; }
    system( $cmd);
}
    
if( ! -s "$namesFilename" ){
    print STDERR "ERROR: NCBI Taxonomy names file \"$namesFilename\" is not found!\n";
    exit(1);
}
if( ! -s "$nodesFilename" ){
    print STDERR "ERROR: NCBI Taxonomy nodes file \"$nodesFilename\" is not found!\n";
    exit(1);
}


# nodes.dmp has the following ranks (with occurrences on the left)
# 273 class
# 6914 family
#  366 forma
# 69791 genus
#   13 infraclass
#   60 infraorder
#    3 kingdom
# 25813 no rank
# 1203 order
#    9 parvorder
#  137 phylum
# 106406 species
#  196 species group
#   81 species subgroup
#  105 subclass
# 1653 subfamily
#  602 subgenus
#    1 subkingdom
#  258 suborder
#   22 subphylum
# 10786 subspecies
#  162 subtribe
#    5 superclass
#  567 superfamily
#    5 superkingdom
#   41 superorder
#    2 superphylum
#  963 tribe
# 4609 varietas

# We want these ranks to appear first in the output.  rank2Index will be extended with all of the ranks (see above) in the file. 
my %rank2Index = ( superkingdom => 0,
		     kingdom => 1,
		     phylum => 2,
		     class => 3,
		     order => 4,
		     family => 5,
		     genus => 6,
		     species => 7 );

my $uniprotSpeciesList = [[]];  # Data structure that mirrors species list file.  Columns are Code	Kingdom	Taxon_Node	Official (scientific) name	Common name	Synonym
my %uniprotSpeciesCodeIndices;  # keys: Uniprot species code; values: index into uniprotSpeciesList
my %speciesCode2taxon; 
my @speciesListHeaders;

my %accession2taxon;

DEBUG( "Reading in Uniprot species list ...");

parseSpecList( $speciesListFilename,    \%speciesCode2taxon, \%uniprotSpeciesCodeIndices, $uniprotSpeciesList, \@speciesListHeaders);
if( $speciesListOldFilename ne "" &&
    -s "$speciesListOldFilename"){
    parseSpecList( $speciesListOldFilename, \%speciesCode2taxon, \%uniprotSpeciesCodeIndices, $uniprotSpeciesList, \@speciesListHeaders);
}

DEBUG( "Reading in uniprot to NCBI Taxon ID mappings from $accession2taxonFilename and putting them into a hash table ...");
open( ACCESSION_2_TAXON, "$accession2taxonFilename") or die( "ERROR: Unable to open accession2taxonFilename \"$accession2taxonFilename\": $!");
my $line;
while( defined( $line = <ACCESSION_2_TAXON>)){
    if( $line !~ m/$accession2taxonRegex/){
	die( "ERROR: Mal-formed line in $accession2taxonFilename: $line (regex: $accession2taxonRegex)");
    }
    
    $accession2taxon{ $1 } = $2;
}
close( ACCESSION_2_TAXON);
# my $uniprotId2TaxonContents = fileToString( $accession2taxonFilename);
# my @uniprotId2TaxonContents = split( /\n/, $uniprotId2TaxonContents);
# for( my $i = 0; $i <= $#uniprotId2TaxonContents; $i++){
#     if( $uniprotId2TaxonContents[ $i] !~ m/^(\S+)\s+(\S+)$/){
# 	die( "ERROR: Mal-formed line in $accession2taxonFilename: $uniprotId2TaxonContents[ $i]");
#     }
#     
#     $accession2taxon{ $1 } = $2;
# }

DEBUG( "Found ". scalar( keys( %accession2taxon)) . " accession to taxon mappings in \"$accession2taxonFilename\"");


DEBUG( "Putting names.dmp into an array, indexed by taxa id...");

# $ head names.dmp
# 1	|	all	|		|	synonym	|
# 1	|	root	|		|	scientific name	|
# 2	|	Bacteria	|	Bacteria <prokaryotes>	|	scientific name	|
# 2	|	Monera	|	Monera <Bacteria>	|	in-part	|


my $namesContents = fileToString( $namesFilename);
my @namesContents = split( /\n/, $namesContents);
my @names; # indexed by taxon id
$#names = 2000000; # reserve some space
for( my $i = 0; $i <= $#namesContents; $i++){
    # Limit search space to scientific names
    if( $namesContents[ $i] =~ m/^(\d+)\s\|\s([^\|]+)\s\|\s.*scientific name/){
	my $taxonId = $1;
	my $name = $2;
	if( $taxonId > $#names){
	    # reserve more memory for the array
	    $#names = $taxonId * 2;
	}
	$names[ $taxonId] = $name;
    }
}
DEBUG( "Largest name index: $#names");


#
DEBUG( "Putting nodes.dmp into two parallel arrays (parent and rank arrays), indexed by taxa id...");
#
# $ head nodes.dmp
# 1	|	1	|	no rank	|		|	8	|	0	|	1	|	0	|	0	|	0	|	0	|	0	|		|
# 2	|	131567	|	superkingdom	|		|	0	|	0	|	11	|	0	|	0	|	0	|	0	|	0	|		|
# 6	|	335928	|	genus	|		|	0	|	1	|	11	|	1	|	0	|	1	|	0	|	0	|		|
# 7	|	6	|	species	|	AC	|	0	|	1	|	11	|	1	|	0	|	1	|	1	|	0	|		|
# 9	|	32199	|	species	|	BA	|	0	|	1	|	11	|	1	|	0	|	1	|	1	|	0	|		|
# 10	|	1706371	|	genus	|		|	0	|	1	|	11	|	1	|	0	|	1	|	0	|	0	|		|

# Columns of interest
# 1,1,no rank
# 2,131567,superkingdom
# 6,335928,genus
# 7,6,species
# 9,32199,species
# 10,1706371,genus
# 11,1707,species
# 13,203488,genus
# 14,13,species
# 16,32011,genus

my $nodesContents = fileToString( $nodesFilename);
my @nodesContents = split( /\n/, $nodesContents);
my @parents; # indexed by taxon id
my @ranks; # indexed by taxon id
$#parents = $#names; # reserve some space
$#ranks   = $#names; # reserve some space

for( my $i = 0; $i <= $#nodesContents; $i++){
    if( $nodesContents[ $i] !~ m/^(\d+)\s+\|\s+(\d+)\s+\|\s+([^\|]+)\s+/){
	die( "ERROR: Mal-formed line (index $i) in $nodesFilename: $nodesContents[ $i]");
    }
    
    my $taxonId = $1;
    my $parent = $2;
    my $rank = $3;
    if( $taxonId > $#parents){
	# reserve more memory for the array
	$#parents = $taxonId * 2;
	$#ranks = $taxonId * 2;
    }
    $parents[ $taxonId] = $parent;
    $ranks[ $taxonId] = $rank;
    if( ! exists( $rank2Index{ $rank})){
	DEBUG( "Found new rank, \"$rank\".  Setting to index " . scalar( keys( %rank2Index)));
	$rank2Index{ $rank} = scalar( keys( %rank2Index));
    }
}
DEBUG( "Largest name index: $#parents");
DEBUG( "Found " . scalar( keys( %rank2Index)) . " unique rank types");

# If an output filename is specified, open it, otherwise, map STDOUT to the OUT file stream
if( $outputFilename ne ""){
    open( OUT, ">$outputFilename") or die("ERROR: Can not open / create output file name \"$outputFilename\": $! ");
}else{
    *OUT = *STDOUT;
}

my @emptyTaxonomy;
print OUT join( $outputSeparator, @speciesListHeaders);  # repeat species list headers
print OUT $outputSeparator . "Accession";

foreach my $taxKey (sort {$rank2Index{$a} <=> $rank2Index{$b}} (keys %rank2Index)){
    print OUT "$outputSeparator$taxKey";
    push( @emptyTaxonomy, "");
}
print OUT "\n";

if( ! $usingSimplifiedAccession2taxon){
    open( ACCESSION_2_TAXON_CACHE, ">$ACCESSION_2_TAXON_ID_FILENAME_DEFAULT") or die( "ERROR: Can not create accession to taxon_id cached file \"$ACCESSION_2_TAXON_ID_FILENAME_DEFAULT\": $!");
}

my $uniprotTuplesContents = fileToString( $uniprotTuplesFilename);
my @uniprotTuplesContents = split( /\n/, $uniprotTuplesContents);
for( my $i = 0; $i <= $#uniprotTuplesContents; $i++){
    my ($accession, $uniprotSpeciesCode) = split( /\t/, $uniprotTuplesContents[ $i]);
    if( $DEBUG >= 2){ DEBUG( "accession: $accession, uniprotSpeciesCode: $uniprotSpeciesCode"); }
    my $taxonId;
    my $key = $uniprotSpeciesCode;  # for caching accession2taxon entries
    if( defined( $uniprotSpeciesCode) &&
	exists( $speciesCode2taxon{ $uniprotSpeciesCode})){
	print OUT join( "$outputSeparator", @{$uniprotSpeciesList->[ $uniprotSpeciesCodeIndices{ $uniprotSpeciesCode} ]}) . "$outputSeparator$accession";
	$taxonId = $speciesCode2taxon{ $uniprotSpeciesCode};
    }elsif( exists( $accession2taxon{ $accession})){
	DEBUG( "NOTE: Found accession to taxon entry for $accession");
	$taxonId = $accession2taxon{ $accession};
	print OUT "$uniprotSpeciesCode$outputSeparator$outputSeparator$taxonId" . "$outputSeparator" x ($#speciesListHeaders + 1 - 3) . "$outputSeparator$accession";
	$key = $accession;
    }else{
	print STDERR "ALERT: Taxon ID not found for $accession in $accession2taxonFilename nor $uniprotSpeciesCode in $speciesListFilename!\n";
	next;
    }

    if( ! $usingSimplifiedAccession2taxon){
	print ACCESSION_2_TAXON_CACHE "$key$outputSeparator$taxonId\n";
    }
    
    my $success = 0;
    my @taxonomy = taxonId2Taxonomy( $taxonId, $uniprotSpeciesCode, $accession, \$success);
    if( $success != 1 &&
	exists( $speciesCode2taxon{ $uniprotSpeciesCode}) &&
	exists( $accession2taxon{ $accession}) &&
	$speciesCode2taxon{ $uniprotSpeciesCode} != $accession2taxon{ $accession}){
	my $taxonIdAccession = $accession2taxon{ $accession};
	DEBUG( "ALERT: Failed to find taxonomy with species code taxon id $taxonId (species code: $uniprotSpeciesCode).  Now trying accession taxon id $taxonIdAccession (accession: $accession) . . .");
	my @taxonomyAccession = taxonId2Taxonomy( $taxonIdAccession, $uniprotSpeciesCode, $accession, \$success);
	if( $success != 1 ){
	    DEBUG( "ALERT: Failed with accession taxon id as well :(");
	}else{
	    @taxonomy = @taxonomyAccession;
	}
    }

    print OUT "$outputSeparator".join("$outputSeparator",@taxonomy)."\n";
}

if( $outputFilename ne ""){
    close(OUT);
}

if( ! $usingSimplifiedAccession2taxon){
    close( ACCESSION_2_TAXON_CACHE);
}

VERBOSE("Done with $0\n");

exit 0;



sub taxonId2Taxonomy{
    my ($taxonId, $uniprotSpeciesCode, $accession, $successRef) = @_;
    if( !defined( $taxonId)){
	die( "ERROR: taxonId2Taxonomy( _BLANK_ )");
    }
    my $currId = $taxonId;
    my @taxonomy = @emptyTaxonomy;

    my $counter = 0;
    $$successRef = ($currId != 1);
    while( $currId != 1){
	if( $counter > 100){
	    die( "ERROR: Too many loops in $nodesFilename (currId: $currId; uniprotSpeciesCode: $uniprotSpeciesCode; accession: $accession; taxonId: $taxonId)");
	}
	if( ! defined( $names[ $currId])){
	    print STDERR "ALERT: No name found for taxa id $currId (uniprotSpeciesCode: $uniprotSpeciesCode; accession: $accession; taxonId: $taxonId)\n";
	    $$successRef = 0;
	    last;
	}
	if( ! defined( $ranks[ $currId])){
	    print STDERR "ALERT: No rank found for taxa id $currId (uniprotSpeciesCode: $uniprotSpeciesCode; accession: $accession; taxonId: $taxonId)\n";
	    $$successRef = 0;
	    last;
	}
	my $rank = $ranks[ $currId];
	$taxonomy[ $rank2Index{ $rank} ] = $names[ $currId];
	#DEBUG( "currId: $currId; name: $names[$currId]; parent: $parents[$currId]; rank: $rank; counter: $counter");
	$currId = $parents[ $currId];
	$counter++;
    }
    return @taxonomy;
}



sub parseSpecList{
    my ( $speciesListFilename, $speciesCode2taxonRef, $uniprotSpeciesCodeIndicesRef, $uniprotSpeciesListRef, $speciesListHeadersRef) = @_;

    # foreach my $header ("Species_Code", "Kingdom", "Taxon_Node", "Official (scientific) name", "Common name", "Synonym"){
    # 	push( @{$speciesListHeadersRef}, $header);
    # }
    if( $#{$speciesListHeadersRef} < 0){
	@{$speciesListHeadersRef} = ( "Species_Code", "Kingdom", "Taxon_Node", "Official (scientific) name", "Common name", "Synonym");
    }
    
    my $uniprotSpeciesCode2TaxonContents = fileToString( $speciesListFilename);
    # remove header and footer
    $uniprotSpeciesCode2TaxonContents =~ s/^.+___+\s+//s;
    $uniprotSpeciesCode2TaxonContents =~ s/\s---+.+$//s;

    # put each common name and/or synonym onto the same line as the official name
    $uniprotSpeciesCode2TaxonContents =~ s/\n\s+(C=|S=)/\t$1/g;
    # Add space for missing common names
    $uniprotSpeciesCode2TaxonContents =~ s/(N=[^\t]+\t)(S=)/$1\t$2/g;
    # remove the (N|C|S)=
    $uniprotSpeciesCode2TaxonContents =~ s/(N|C|S)=//g;

    # finish changing to tab-delimited (and ignore the ":")
    $uniprotSpeciesCode2TaxonContents =~ s/^(\S+)\s+(\w)\s+(\d+):?\s+/$1\t$2\t$3\t/gm;


    # remove "virtual codes" section break
    # =======================================================================
    # (2) "Virtual" codes that regroup organisms at a certain taxonomic level
    # =======================================================================
    $uniprotSpeciesCode2TaxonContents =~ s/\s*===+\s+[^\n]+\s+===+\s*/\n/s;
    
    my @uniprotSpeciesCode2TaxonContents = split( /\n/, $uniprotSpeciesCode2TaxonContents);
    DEBUG("Found " .($#uniprotSpeciesCode2TaxonContents+1)." entries in $speciesListFilename");
    my $codeIndex = $#{$uniprotSpeciesListRef} + 1;  # start where we left off if a specList has already been parsed or at 0 if it's the first time
    for( my $i = 0; $i <= $#uniprotSpeciesCode2TaxonContents; $i++){
	my @line = split( /\t/, $uniprotSpeciesCode2TaxonContents[ $i]);
	if( ! exists( $uniprotSpeciesCodeIndices{ $line[0] })){
	    $uniprotSpeciesCodeIndices{ $line[0] } = $codeIndex;  # Using the species code as the key, save the index of this entry for faster look-up
	    $speciesCode2taxon{ $line[$SPECIES_LIST_COL_CODE] } = $line[ $SPECIES_LIST_COL_TAXON];  # Using the species code as the key, save the index of this entry for faster look-up
	    my $j = 0;
	    for( $j = 0; $j <= $#line; $j++){
		$uniprotSpeciesListRef->[ $codeIndex][$j] = $line[$j];
	    }
	    # add blank entries if necessary
	    for( ; $j <= $#{$speciesListHeadersRef}; $j++){
		$uniprotSpeciesListRef->[ $codeIndex][$j] = "";
	    }

	    $codeIndex++;
	}else{
	    # DEBUG( "Note: Already found a speclist entry for $line[0] (while parsing $speciesListFilename)"); 
	}
    }

    DEBUG( "Now there are $codeIndex rows in the speclist data structure (after reading in species list file \"$speciesListFilename\")");
}


sub downloadFile{
    my ($url, $filename) = @_;
    
    my $which_curl = `which curl`;
    my $which_wget = `which wget`;
    if( $which_curl ne ""){
	DEBUG( "Downloading $url with curl");
	my $rc = system( "curl $url -o $filename");
	if( $rc != 0){
	    die("ERROR: Failed retrieving $url (with curl)");
	}
    }elsif( $which_wget ne ""){
	DEBUG( "Downloading $url with wget");
	my $rc = system( "wget $url -O $filename");
	if( $rc != 0){
	    die("ERROR: Failed retrieving $url (with wget)");
	}
    }else{
	die("ERROR: Could not find curl nor wget to download $filename!");
    }

    # test that the file now exists
    if( ! -s "$filename" ){
	die( "ERROR: Attempt to download \"$filename\" from \"$url\" failed!");
    }
}


__END__

=head1 NAME

uniprot2taxonomy.pl - For each specified uniprot id, get the taxon id and its taxonomy

=head1 SYNOPSIS

uniprot2taxonomy.pl INPUT [OPTIONS]

 Input: 
   --uniprots <filename with accessions or tab-delimited accessions and Uniprot species codes>
    
 Options (all are optional):
   --tempDir <temporary directory name for cache files (defaults to current directory)>
   --speclist <Uniprot species codes information filename (input or the filename for a cache file) (defaults to tempDir/speclist.tab)>
   --speclistOld <Uniprot species codes information filename>
   --names <NCBI Taxonomy names filename (input or the filename for a cache file) (defaults to tempDir/names.dmp)>
   --nodes <NCBI Taxonomy nodes filename (input or the filename for a cache file) (defaults to tempDir/nodes.dmp)>
   --uniprotMappings <accession to taxon mappings filename (input or the filename for a cache file) (defaults to tempDir/idmapping_selected-cached.tab)>
   --out <taxonomy output filename> (defaults to STDOUT)
   --separator <separator> (default to tab character)
   -v | --verbose  Display verbose output
   -d | --debug  Display debugging information
   --help  Brief help message
   --man  Full documentation
=cut
