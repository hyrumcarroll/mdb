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
# Copyright 2015-2017 Hyrum D. Carroll

# The term "Index" in a variable name identifies the location (index value: 0 .... len - 1) of the start, end, etc. of a string or array
# Iterates through one or more spouge format retrievals and truncates list with a score (E-value) greater than the user specified threshold

use strict;
use warnings;

my $FALSE = 0;
my $TRUE = 1;

#my $E_VALUE_THRESHOLD = 0.000001;


my $usage = "Usage: $0:
  -spouge=<SPOUGE output file>
  -spougee=<SPOUGEe output file>
  -spougeeEValue=<Threshold E-Value for SPOUGEe>
  [-v]
  [-d]";

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
my $spougeOutputFileName = "";
my $spougeeOutputFileName = "";
my $spougeeEValue = "";

if( exists($args{"v"}) ||
    exists($args{"verbose"})){
    $VERBOSE = $TRUE;
}

if( exists($args{"d"}) ||
    exists($args{"debug"})){
    $DEBUG = $TRUE;
    $VERBOSE = $TRUE;
}

if( defined( $args{"spouge"})){
    $spougeOutputFileName = $args{"spouge"};
}
if( $spougeOutputFileName eq ""){
    print STDERR "ERROR: SPOUGE output file not specified!\n";
    die( "$usage\n");
}

if( defined( $args{"spougee"})){
    $spougeeOutputFileName = $args{"spougee"};
}
if( $spougeeOutputFileName eq ""){
    print STDERR "ERROR: SPOUGEe output file not specified!\n";
    die( "$usage\n");
}

if( defined( $args{"spougeeevalue"})){
    $spougeeEValue = $args{"spougeeevalue"};
}
if( $spougeeEValue eq ""){
    print STDERR "ERROR: SPOUGEe E-Value not specified!\n";
    die( "$usage\n");
}

#open(SPOUGE, "$spougeOutputFileName") or die( "ERROR: Can not open / create SPOUGE output file \"$spougeOutputFileName\": $! ");

my $spougeFile = fileToString($spougeOutputFileName);
#print STDERR "length(\$spougeFile): " . length($spougeFile) . "\n";
# split up spouge file into individual retrievals (based on 2+ \n's)
# foreach individual retrieval
# + get header
# + get total # of relevant records
# + include all records <= threshold
my @retrievals = split( /\n\n+/, $spougeFile);
#print STDERR "\$#retrievals: $#retrievals\n";
open(SPOUGEE, ">$spougeeOutputFileName") or die( "ERROR: Can not open / create SPOUGEe output file \"$spougeeOutputFileName\": $! ");
for( my $retrievalI = 0; $retrievalI <= $#retrievals; $retrievalI++){
    #print STDERR "retrievalI: $retrievalI\n";
    my @lines = split( /\n/, $retrievals[$retrievalI]);
    #print STDERR "\t\$#lines: $#lines\n";
    
    # header
    my $header = $lines[0];
    print SPOUGEE "$header\n";

    $header = $lines[1];
    if( $header !~ m/^(\d+)\s*/){
	die("ERROR: Mal-formed header \"$header\".  Expecting total number of relevant records in the 1st column ");
    }
    print SPOUGEE "$header\n";

    for( my $lineI = 2; $lineI <= $#lines; $lineI++){
	if( $lines[$lineI] =~ m/^\s+$/){
	    last;
	}
	if( $lines[$lineI] !~ m/^([01]\s+)([\d\.\+\-eE]+)/){
	    chomp;
	    die( "ERROR: Mal-formed line: \"$lines[$lineI]\".  Expecting E-Value in the 2nd column ");
	}
	if( $2 > $spougeeEValue){
	    last;
	}
	print SPOUGEE "$lines[$lineI]\n";
    }
    print SPOUGEE "\n";
}
#close( SPOUGE);
close( SPOUGEE);

exit(0);


sub fileToString {
            my $file = $_[0];
            
            my $terminator = $/;
            undef $/;
            open INPUT, "<$file" or die("Input file '" . $file . "' failed to open.\n" . "Died");
            my $str = <INPUT>; # terminator undefined : $str is the whole file.
            close INPUT;
        
            $/ = $terminator; # Restore for normal behavior later
        
            return $str;
        }
