#!/usr/bin/perl

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

BEGIN
{
    # Look for perl modules in ./perl
    unshift( @INC, "./perl"); 
}

use HDC::Common qw(:ALL);


if( $#ARGV + 1 < 2){
    die( "Usage $0:  <randomQueriesFile (e.g., queriesPerDA-diverse.tab)>  <multiDomainQueriesFile (e.g., domainQueries_multi.tab)>");
}

my $randomQueriesFile=$ARGV[0];
my $multiDomainQueriesFile=$ARGV[1];

my @multiDomainQueriesFileContents = split( "\n", fileToString($multiDomainQueriesFile));
my %taxonInfo;

for( my $line = 1; $line <= $#multiDomainQueriesFileContents; $line++){  # ignore the header
    if( $multiDomainQueriesFileContents[$line] !~ m/^(\S+)\s+(.+)$/){
	die( "Mal-formed domain queries info file line: \"$multiDomainQueriesFileContents[$line]\"!");
    }
    my $taxon = $1;
    my $info = $2;
    
    $taxonInfo{$taxon} = $info;
}

my @randomQueriesFileContents = split("\n", fileToString( $randomQueriesFile));

# remove header line if it's a comment
if( $randomQueriesFileContents[0] =~ m/^\s*#/){
    shift @randomQueriesFileContents;
}

# get random taxon per each DA
# DA_ID	Random_Taxon	Shortest_Taxon	Longest_Taxon
my @randomTaxon;
foreach my $line (@randomQueriesFileContents){
    my @line = split( "\t", $line);
    push( @randomTaxon, $line[1]);
}

print "$multiDomainQueriesFileContents[0]\n";
foreach my $taxon (@randomTaxon){
    if( ! exists( $taxonInfo{$taxon})){
	die( "ERROR: Could not find taxon \"$taxon\" in $multiDomainQueriesFile!");
    }
    print "$taxon\t$taxonInfo{$taxon}\n";
}

exit
