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

# Description: For each of the domain architecture (DA) specied in the first
#              file (e.g., a list a queries), count the number of sequences that
#              are members of that same DA.

use strict;
use warnings;

BEGIN
{
    # Look for perl modules in ./perl
    unshift( @INC, "./perl"); 
}

use HDC::Common qw(:ALL);

if( $#ARGV + 1 < 2){
    die( "Usage $0:  <queries filename (e.g., queriesPerDA-diverse.tab)>  <relevance file (e.g., taxon2da.tab)>");
}

my $queriesFile=$ARGV[0];
my $relevanceFile=$ARGV[1];

my @relevanceFileContents = split( "\n", fileToString($relevanceFile));
my %daCounts;
for( my $line = 1; $line <= $#relevanceFileContents; $line++){  # ignore the header
    if( $relevanceFileContents[$line] !~ m/^(\S+)\s+(\S+)/){
	die( "Mal-formed relevance info file line: \"$relevanceFileContents[$line]\"!");
    }
    my $taxon = $1;
    my @das = split( ',', $2);

    for( my $daI = 0; $daI <= $#das; $daI++){
	if( ! exists($daCounts{$das[$daI]})){
	    $daCounts{$das[$daI]} = 1;
	}else{
	    $daCounts{$das[$daI]}++;
	}
    }
}
DEBUG( "Found " . (keys %daCounts) . " DA entries in $relevanceFile");

my @queriesFileContents = split("\n", fileToString( $queriesFile));

DEBUG( "Found " . ($#queriesFileContents + 1) . " DA entries in $queriesFile");

# remove header line if it's a comment
if( $queriesFileContents[0] =~ m/^\s*#/){
    shift @queriesFileContents;
}

# get random taxon per each DA
# DA_ID	Random_Taxon	Shortest_Taxon	Longest_Taxon
my @das_subset;
foreach my $line (@queriesFileContents){
    my @line = split( "\t", $line);
    push( @das_subset, $line[0]);
}

print "#DA,DA_size\n";
foreach my $da (@das_subset){
    if( ! exists( $daCounts{$da})){
	die( "ERROR: Could not find domain architecture \"$da\" in $relevanceFile!");
    }
    print "$da,$daCounts{$da}\n";
}

exit
