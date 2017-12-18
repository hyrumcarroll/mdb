#!/usr/bin/env perl

#############################################################################
# Description: Remove taxa that match the specified labels
# History:     Spawned from chooseRandomSequences.pl
#############################################################################

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

use Getopt::Long;

BEGIN
{
    # Look for perl modules in ./perl
    unshift( @INC, "./perl"); 
}

use HDC::Common qw(:ALL);
use HDC::hyBioScripts;

#
# Command line options
#
my $testOnly = 0;
my $filenameIn;
my $filenameOut;  # = ""
#my $outputType = "";
my $remove;
my $removeFilename;
my $seed;
my $man = 0;
my $help = 0;


my $usageStr = "Usage: $0\n".
"-in <input file name>\n".
"-remove <LABEL> | -removeFile <filename of label(s) to remove>\n".
"[-out <output file name> (defaults to STDOUT)]\n".
"-verbose\n".
"-debug\n".
"-help";


unless (GetOptions("verbose+"            => \$VERBOSE,
		   "debug+"              => \$DEBUG,
		   "in=s"                => \$filenameIn,
		   "out=s"               => \$filenameOut,
		   "remove=s"            => \$remove,
		   "removeFile=s"        => \$removeFilename,
		   'help|?'              => \$help,
		   'man'                 => \$man)){
    print STDERR "$usageStr\n";
    exit (0);
}

if( ! defined( $filenameIn) ||
    (! defined( $remove) &&
     ! defined( $removeFilename))){
    die( "$usageStr\n");
}

my @taxonNames;
my @seqs;

HDC::hyBioScripts::readDataFile( $filenameIn, \@taxonNames, \@seqs);

my @labelsToRemove;

if( defined( $remove)){
    push( @labelsToRemove, $remove);
}else{
    my $labelsFileContents = fileToString( $removeFilename);
    @labelsToRemove = split( /\s+/, $labelsFileContents);
    DEBUG( "\@labelsToRemove: @labelsToRemove");
}

DEBUG( ($#taxonNames + 1) . " taxa before removing any");
HDC::hyBioScripts::removeTaxa(\@taxonNames, \@seqs, \@labelsToRemove);
DEBUG( ($#taxonNames + 1) . " taxa after removing ". ($#labelsToRemove + 1)." taxa");

HDC::hyBioScripts::printDataSet(\@taxonNames, \@seqs, $filenameOut); #, $type);

exit(0);
