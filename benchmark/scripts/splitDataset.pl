#!/usr/bin/env perl

#############################################################################
# Description: Split each sequence in a specified dataset into it's own file.
#              If a file of taxa labels is supplied, only include those
#              sequences.
#              Equivalent to:
#              for taxon in `cat $labelsFilename`; do
#                  taxonFilename="${taxon//|/_}"  # replace all occurrences of "|" with "_"
#                  fastaFilename="$directoryOut/$taxonFilename.fa"
#                  grep -A 1 "^>$taxon"  "$filenameIn" > "$fastaFilename"
#              done
# History:     Spawned from removeTaxa.pl
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
# Copyright 2016-2017 Hyrum D. Carroll

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
my $directoryOut = "./";  #default to current directory
my $labelsFilename;
#my $outputType = "";
my $seed;
my $man = 0;
my $help = 0;


my $usageStr = "Usage: $0\n".
"-in <input file name>\n".
#"[-out <output file name> (defaults to STDOUT)]\n".
"[-labels <filename of label(s)>\n".
"-verbose\n".
"-debug\n".
"-help";


unless (GetOptions("verbose+"            => \$VERBOSE,
		   "debug+"              => \$DEBUG,
		   "in=s"                => \$filenameIn,
		   "dir=s"               => \$directoryOut,
		   "labels|label=s"      => \$labelsFilename,
		   'help|?'              => \$help,
		   'man'                 => \$man)){
    print STDERR "$usageStr\n";
    exit (0);
}

if( ! defined( $filenameIn)){
    die( "$usageStr\n");
}

my @taxonNames;
my @seqs;

HDC::hyBioScripts::readDataFile( $filenameIn, \@taxonNames, \@seqs);

my %labels;  # if populated, only process these labels

if( defined( $labelsFilename)){
    my $labelsFileContents = fileToString( $labelsFilename);
    my @labels = split( /\s+/, $labelsFileContents);
    DEBUG( "\@labels: @labels");
    for my $label (@labels){
	$labels{$label} = 1;
    }
}

# foreach taxon in the dataset, check if it's in the labels file (if it exists)
for( my $i = 0; $i <= $#taxonNames; $i++){

    # if the labels hash is empty or we found the taxon
    $taxonNames[$i] =~ s/\|/_/g;  # replace all occurrences of "|" with "_"
    if( ! %labels || exists( $labels{ $taxonNames[$i] })){
	my @taxonName = ( $taxonNames[$i]);
	my @seq = ( $seqs[$i] );
	my $filename = "$directoryOut/$taxonNames[$i].fa";
	HDC::hyBioScripts::printDataSet(\@taxonName, \@seq, $filename); #, $type);
    }
}


exit(0);
