#!/usr/bin/perl

################################################################################
# Description: See POD documentation at the end of the file
################################################################################

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

#
# Command line options
#
# Configuration defaults, can be overriden using command line arguments
my $testOnly = 0;
my $mainFilenameIn = "";
my $keysFilenameIn = "";
my $filenameOut = "";
my $separator = "\t"; # default to tab-delimited
my $sorted = $FALSE;
my $man = 0;
my $help = 0;


unless (GetOptions("verbose+"            => \$VERBOSE,
		   "debug+"              => \$DEBUG,
		   "in|main=s"           => \$mainFilenameIn,
		   "keys=s"              => \$keysFilenameIn,
		   "out=s"               => \$filenameOut,
		   "sorted=s"            => \$sorted,
		   "separator=s"         => \$separator,
		   "testOnly"            => \$testOnly,
		   'help|?'              => \$help,
		   'man'                 => \$man)){
    #print STDERR "$usageStr\n";
    #exit (0);
    pod2usage(2);
}

# if -in, -keys and or -out are missing, look at ARGV
if( ! defined( $mainFilenameIn) ||
    $mainFilenameIn eq ""){
    if( $#ARGV >= 0){
	$mainFilenameIn = shift @ARGV;
    }else{
	print STDERR "\nERROR: Main input file required (use -in)\n\n";
	pod2usage(2);
    }
}

if( ! defined( $keysFilenameIn) ||
    $keysFilenameIn eq ""){
    if( $#ARGV >= 0){
	$keysFilenameIn = shift @ARGV;
    }else{
	print STDERR "\nERROR: Keys input file required (use -keys)\n\n";
	pod2usage(2);
    }
}

if( ! defined( $filenameOut) ||
    $filenameOut eq ""){
    if( $#ARGV >= 0){
	$filenameOut = shift @ARGV;
    }
}

DEBUG( "$0 parameters:
    mainFilenameIn:     $mainFilenameIn
    keysFilenameIn:     $keysFilenameIn
    filenameOut:        $filenameOut
    sorted:             $sorted
    separator:          $separator
    testOnly:           $testOnly
    help:               $help
    man:                $man
    VERBOSE:            $VERBOSE
    DEBUG:              $DEBUG
");

pod2usage(1) if $help;
pod2usage(-exitval => 0, -verbose => 2) if $man;

#
# OUTLINE
#

# open main file
# open keys file
# open output file

# read in 1 line from main file and key file
# echo main file header(s)
# ignore keys file header(s)

# while main file is not EOF and keys file is not EOF:
#   compare the first columns
#   if the main file is  < key, advance the main file
#   if the main file is  > key, advance the keys file
#   if the main file is == key, echo the main file and advance both files

# close open files


# open main file
open( MAIN, "$mainFilenameIn") or die "ERROR: Can not open (main) input file \"$mainFilenameIn\": $! ";

# open keys file
open( KEYS, "$keysFilenameIn") or die "ERROR: Can not open (keys) input file \"$keysFilenameIn\": $! ";

# read in 1 line from main file and key file

# echo main file header(s)
my $mainFileLine;
my $mainHeaders = "";
# echo headers, if any (but ONLY the first contiguous lines starting with #)
$mainFileLine = <MAIN>;
while( defined( $mainFileLine) && $mainFileLine =~ m/^\s*#/){
    $mainHeaders .= $mainFileLine;
    chomp $mainFileLine;
    DEBUG("Found header mainFileLine \"$mainFileLine\" from $mainFilenameIn.");
    $mainFileLine = <MAIN>;
}

# ignore keys file header(s)
my $keysFileLine;
# echo headers, if any (but ONLY the first contiguous lines starting with #)
$keysFileLine = <KEYS>;
while( defined( $keysFileLine) && $keysFileLine =~ m/^\s*#/){
    chomp $keysFileLine;
    DEBUG("Found header keysFileLine \"$keysFileLine\" from $keysFilenameIn.");
    $keysFileLine = <KEYS>;
}

if( ! defined( $mainFileLine) ||
    ! defined( $keysFileLine)){
    close( MAIN);
    close( KEYS);
    die( "ALERT: Either the main file and/or the keys file had no (real) content!");
}

# open output file
if( $filenameOut ne ""){
    open( OUT, ">$filenameOut") or die "ERROR: Can not open output file \"$filenameOut\": $! ";
}else{

    # # make STDOUT "HOT" (unbuffered)
    # my $ofh = select STDOUT;
    # $| = 1;
    # select $ofh;

    DEBUG("Using STDOUT in lieu of output file");
    *OUT = *STDOUT;
}

if( $mainHeaders ne ""){
    print OUT $mainHeaders;
}

if( $sorted){
    # while main file is not EOF and keys file is not EOF:
    while( $TRUE){
	#   compare the first columns
	if( $mainFileLine !~ m/^\s*([^$separator]+)/){
	    chomp $mainFileLine;
	    die("ERROR: Malformed main file line: \"$mainFileLine\"");
	}
	my $mainFirstColumn = $1;
	chomp $mainFirstColumn;
	
	if( $keysFileLine !~ m/^\s*([^$separator]+)/){
	    chomp $keysFileLine;
	    die("ERROR: Malformed keys file line: \"$keysFileLine\"");
	}
	my $keysFirstColumn = $1;
	chomp $keysFirstColumn;

	DEBUG( "mainFirstColumn: $mainFirstColumn; keysFirstColumn: $keysFirstColumn;");

	if( $mainFirstColumn lt $keysFirstColumn){
	    #   if the main file is  < key, advance the main file
	    if( ! defined( $mainFileLine = <MAIN>)){
		last;
	    }

	}elsif( $mainFirstColumn gt $keysFirstColumn){
	    #   if the main file is  > key, advance the keys file
	    if( ! defined( $keysFileLine = <KEYS>)){
		last;
	    }

	}else{
	    #   if the main file is == key, echo the main file and advance both files
	    print OUT $mainFileLine;
	    if( ! defined( $mainFileLine = <MAIN>) ||
		! defined( $keysFileLine = <KEYS>)){
		last;
	    }
	}
    }
}else{

    #
    # unsorted, hash key values
    #
    # hash first column of keys file
    # iterate through the main file, printing out matches

    my %keys;
    do{
	if( $keysFileLine !~ m/^\s*([^$separator]+)/){
	    chomp $keysFileLine;
	    die("ERROR: Malformed keys file line: \"$keysFileLine\"");
	}
	my $keysFirstColumn = $1;
	chomp $keysFirstColumn;
	$keys{$keysFirstColumn} = 1;
    }while( defined( $keysFileLine = <KEYS>));

    DEBUG( "Found " . scalar(keys( %keys)). " entries from $keysFilenameIn");


    # iterate through the main file, printing out matches
    do{
	if( $mainFileLine !~ m/^\s*([^$separator]+)/){
	    chomp $mainFileLine;
	    die("ERROR: Malformed main file line: \"$mainFileLine\"");
	}
	my $mainFirstColumn = $1;
	chomp $mainFirstColumn;
	if( defined( $keys{$mainFirstColumn})){
	    print OUT $mainFileLine;
	# }else{
	#     DEBUG("$mainFirstColumn not found in $keysFilenameIn");
	}
    }while( defined( $mainFileLine = <MAIN>));
}

# close open files
close( MAIN);
close( KEYS);
if( $filenameOut ne ""){
    close( OUT);
}

exit(0);


__END__
    
=head1 NAME

filterByUnion.pl
    
=head1 SYNOPSIS

  -in|-main <main filename>
  -keys <keys filename>
  [-out <output filename>]
  [options]

=head1 OPTIONS

=over 8

=item B<-in|-main> <main input file name>

Filename for file whose lines will be displayed

=item B<-keys> <keys input file name>

Only the first column of this file will be used to check for a union in the main input file

=item B<-out> <output file name>

Defaults to STDOUT

=item B<-sorted>

If both files are sorted, then traverse the files linearly.  Otherwise, hash the key values. Defaults to not set.

=item B<-separator> <separator>

Defaults to $separator

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-v|verbose>

Defaults to $VERBOSE

=item B<-d|debug>

Defaults to $DEBUG

=back

=head1 DESCRIPTION

Output the lines of a file, where the first column matches the first column in a second file.

Note: Assumes the two files are sorted (for performance reasons).
    
Note: Comparison works for strings and equal length numerical values.

=cut
