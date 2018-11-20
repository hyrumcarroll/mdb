#!/usr/bin/perl

use strict;
use warnings;

#
# Overview
#

# Get domain locations filename
# Get FASTA filename

# Put each of the FASTA sequences into a dictionary with the label as the key

# Create a single-domain FASTA file
# Create a relevancy file

# Parse the domain locations file
#     Split the file by ">"
#         Split each of the records by "\n"
#         Get the label from the first line
#         Get the corresponding sequence
#         For each of the rest of the lines
#             Get domain name and start and end of the domain
#             Create a new label, based on the record's label and the domain index
#             Create a new (single-domain) sequence
#             Add an entry in the relevance database: new label and the domain name


sub fileToString{
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


# Get domain locations filename
# Get FASTA filename

my $domainLocsFilename = shift || "../benchmark/domainLocs.tab";
my $fastaInFilename = shift || "../benchmark/library_all_domains_MDB.fa";

my $fastaOutFilename = shift || "singleDomains.fa";
my $relevencyFilename = shift || "relevancy.tab";
    
print "domainLocsFilename: $domainLocsFilename\n";
print "fastaInFilename: $fastaInFilename\n";

# Put each of the FASTA sequences into a dictionary with the label as the key
my %fasta;
my @fastaRecords = split( />/, fileToString( $fastaInFilename));
my $firstRecord = shift @fastaRecords;
print "" . ($#fastaRecords + 1) . " FASTA records found (ignoring first record: $firstRecord)\n";
foreach my $fasta (@fastaRecords){
    #print "fasta: $fasta\n";
    my @fastaFields = split( /\n/, $fasta);
    if( $#fastaFields < 1){
	die("fastaFields: @fastaFields");
    }
    $fasta{ $fastaFields[0] } = join('', @fastaFields[1..$#fastaFields]);
    #print STDERR "Added $fastaFields[0]: " . join('', @fastaFields[1..$#fastaFields]) ."\n";
}
print( "".scalar (keys %fasta)." FASTA records processed\n");



# Create a single-domain FASTA file
# Create a relevancy file

open( FASTA, ">$fastaOutFilename") or die("ERROR: Unable to create fastaOutFilename \"$fastaOutFilename\": $!");
open( HOMOLOGY, ">$relevencyFilename") or die("ERROR: Unable to create relevencyFilename \"$relevencyFilename\": $!");
print HOMOLOGY "# Label\tdomain_name\n";


# >up_Q5LL22_Q5LL22_SILPO
# CL3	2	66
# CL276	76	262
# CL20	661	694
# >up_Q5LL78_Q5LL78_SILPO
# CL123	2	65
# CL61	82	431
# >up_Q5LL87_Q5LL87_SILPO


# Parse the domain locations file
# Split the file by ">"
my @domainsRecords = split( />/, fileToString($domainLocsFilename));
# Split each of the records by "\n"
print "Found $#domainsRecords domains records\n";
my $firstDomainsRecord = shift @domainsRecords;  # should be empty
if( length($firstDomainsRecord) > 0){
    die("firstDomainsRecord: $firstDomainsRecord");
}
foreach my $domainsRecord (@domainsRecords){
    my @domainsFields = split(/\n/, $domainsRecord);
    # Get the label from the first line
    my $domainsLabel = shift @domainsFields;
    # Get the corresponding sequence
    if( ! defined( $fasta{ $domainsLabel })){
	die("domainsLabel: $domainsLabel\n");
    }
    my $seq = $fasta{ $domainsLabel};

    
    # print "domainsFields: @domainsFields\n";
    
    # For each of the rest of the lines
    my $index = 0;
    my @domainNames;
    my @startPositions;
    my @endingPositions;
    
    foreach my $domainRecord (@domainsFields){
        # Get domain name and start and end of the domain
	($domainNames[$index], $startPositions[$index], $endingPositions[$index]) = split( /\s/, $domainRecord);
	#print "Found $domainName, $startPos, $endingPos\n";
	$index++;
    }

    my $previousEndingPosition = 0;
    for( my $i = 0; $i <= $#domainNames; $i++){
        # Create a new label, based on the record's label and the domain index
        my $newLabel = "${domainsLabel}_$i";

	my $nextStartingPosition = length( $seq) + 1;
	if( $i != $#domainNames ){
	    $nextStartingPosition = $startPositions[$i + 1 ];
	}
	
        # Create a new (single-domain) sequence
	print FASTA ">$newLabel\n" . substr( $seq, $previousEndingPosition, $nextStartingPosition - $previousEndingPosition - 1) . "\n";
	
        # Add an entry in the relevance database: new label and the domain name
	print HOMOLOGY "$newLabel\t$domainNames[$i]\n";

	$previousEndingPosition = $endingPositions[$i];
    }
}

close( HOMOLOGY );
close( FASTA );
