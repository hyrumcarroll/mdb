package HDC::Common;

use strict;
use warnings;

use Exporter;

our $VERSION     = 1.00;
our @ISA         = qw(Exporter);
our @EXPORT      = ();
our @EXPORT_OK   = qw($FALSE $TRUE $VERBOSE $DEBUG fileToString VERBOSE DEBUG MIN MAX); # put stuff here you want to export
our %EXPORT_TAGS = ( DEFAULT => [qw($FALSE $TRUE $VERBOSE $DEBUG)],
		     ALL     => [@EXPORT_OK]);

=head 1

# Functions and globals common to all HDC scripts
=cut

#
# GLOBALS
#
our $FALSE = 0;
our $TRUE = 1;

our $VERBOSE = 0;
our $DEBUG = 0;

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


1;
