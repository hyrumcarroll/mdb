#!/usr/bin/env perl 

################################################################################
# File:         fileSplicer.pl
# Author:       Hyrum Carroll
# Description:  Breaks up a file into $pieces (with each subsets' filename being the original's with a number appended to it)
# Example:      An input file with consecutive numbers 1-15, and 3 pieces:
#               file.0: 1,4,7,10,13
#               file.1: 2,5,8,11,14
#               file.2: 3,6,9,12,15
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

use warnings;
use strict;

my $DEFAULT_NUM_PIECES = 10;


# make STDOUT "HOT" (unbuffered)
my $ofh = select STDOUT;
$| = 1;
select $ofh;

my $FALSE = 0;
my $TRUE = 1;


my $VERBOSE = $FALSE;
my $DEBUG = $FALSE;

my $file = shift;
if( ! defined ($file)){ die( "Usage $0: <file> [ <number of subsets> ]\n"); }
my $pieces = shift;
if( ! defined( $pieces)){
    $pieces = $DEFAULT_NUM_PIECES;
}

if( $VERBOSE || $DEBUG){
    print STDERR "Config:\n";
    print STDERR "   file           $file\n";          
    print STDERR "   pieces         $pieces\n";          
    print STDERR "   VERBOSE        $VERBOSE\n";     
    print STDERR "   DEBUG          $DEBUG\n";          
}

# read in the file
open(FILE, $file) or die( "Error opening file \"$file\": $! ");
my @lines = <FILE>;
close( FILE);

my $numLines = $#lines + 1;

my $subsetFileName = "";
my $linesPerSubset = ($#lines + 1) / $pieces;
if( $VERBOSE || $DEBUG){
    print STDERR "linesPerSubset: $linesPerSubset\n";
}

for( my $i = 0; $i < $pieces; $i++){
    if( $VERBOSE || $DEBUG){
	print STDERR "Subset: $i\n";
    }
    $subsetFileName = "$file.$i";
    open( SUBSET, ">$subsetFileName") or die( "ERROR: Can not open subset file \"$subsetFileName\": $! ");
    
    for( my $lineIndex = $i; $lineIndex < $numLines; $lineIndex += $pieces){
	print SUBSET $lines[$lineIndex];
    }
    close( SUBSET);
}

exit (0);

