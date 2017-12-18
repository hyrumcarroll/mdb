package HDC::SortingShuffling;

use strict;
use warnings;

use HDC::Common qw(:ALL);

use Exporter;

our $VERSION     = 1.00;
our @ISA         = qw(Exporter);
our @EXPORT      = ();
our @EXPORT_OK   = qw(fisherYatesShuffle sort2Arrays merge2arrays mergesort2arrays quickSort2Arrays quickSort2Arrays); # put stuff here you want to export
our %EXPORT_TAGS = ( DEFAULT => [],
		     ALL     => [@EXPORT_OK]);

=head 1

# Sorting and shuffling functions
=cut

# randomly permutate @array in place
sub fisherYatesShuffle{
    my $arrayRef = shift;
    my $i = @{$arrayRef};
    while( --$i ){
        my $j = int(rand( $i+1 ));
        @{$arrayRef}[$i,$j] = @{$arrayRef}[$j,$i];
    }
}



sub sort2Arrays{
    my ($keysRef, $otherRef) = @_;

    if( $#{$keysRef} < $#{$otherRef}){
	die( "$#{$keysRef} < $#{$otherRef}!; ");
    }
    
    my @tempArrayA;
    my @tempArrayB;
    $#tempArrayA = $#{$keysRef};
    $#tempArrayB = $#{$otherRef};

    mergesort2arrays( $keysRef, $otherRef, 0, $#{$keysRef}, \@tempArrayA, \@tempArrayB);
}


sub merge2arrays{
    my ($keysRef, $otherRef, $startIndex, $endIndex, $tempKeysArrayRef, $tempOtherArrayRef) = @_;
    #print STDERR "merge2arrays( keysRef, otherRef, $startIndex, $endIndex, tempKeysArrayRef, tempOtherArrayRef)\n";

    my $midIndex = int( ($startIndex + $endIndex) / 2); # index of midpoint
    my $first1 = $startIndex;       # beginning of first subarray
    my $last1  = $midIndex;         # end of first subarray
    my $first2 = $midIndex + 1;    # beginning of second subarray
    my $last2  = $endIndex;        # end of second subarray

    #
    # while both subarrays are not empty, copy the smaller item into the temporary array
    #
    my $index = $first1;    # next available location in tempArray
    for (; ($first1 <= $last1) && ($first2 <= $last2); ++$index){
	# Invariant: tempArray[first..index-1] is in order
	if( $keysRef->[$first1] <= $keysRef->[$first2]){
	    $tempKeysArrayRef->[$index] = $keysRef->[$first1];
	    $tempOtherArrayRef->[$index] = $otherRef->[$first1];
	    ++$first1;
	}else{
	    $tempKeysArrayRef->[$index] = $keysRef->[$first2];
	    $tempOtherArrayRef->[$index] = $otherRef->[$first2];
	    ++$first2;
	}  # end if
    }  # end for

    #
    # finish off the nonempty subarray
    #
    # finish off the first subarray, if necessary
    for( ; $first1 <= $last1; ++$first1, ++$index){
	# Invariant: tempArray[first..index-1] is in order
	$tempKeysArrayRef->[$index] = $keysRef->[$first1];
	$tempOtherArrayRef->[$index] = $otherRef->[$first1];
    }

    # finish off the second subarray, if necessary
    for( ; $first2 <= $last2; ++$first2, ++$index){
	# Invariant: tempArray[first..index-1] is in order
	$tempKeysArrayRef->[$index] = $keysRef->[$first2];
	$tempOtherArrayRef->[$index] = $otherRef->[$first2];
    }

    #
    # copy the result back into the original array
    #
    for( $index = $startIndex; $index <= $endIndex; ++$index){
	$keysRef->[$index] = $tempKeysArrayRef->[$index];
	$otherRef->[$index] = $tempOtherArrayRef->[$index];
    }

    #print STDERR "\tkeys: @{$keysRef}\n";
    #print STDERR "\tother: @{$otherRef}\n";
}

sub mergesort2arrays{
    my ($keysRef, $otherRef, $startIndex, $endIndex, $tempKeysArrayRef, $tempOtherArrayRef) = @_;
    #print STDERR "mergesort2arrays( keysRef, otherRef, $startIndex, $endIndex, tempKeysArrayRef, tempOtherArrayRef)\n";

    # if( $#{$keysRef} != $#{$otherRef}){
    # die( "$#{$keysRef} != $#{$otherRef}!; ");
    # }

    if( $endIndex - $startIndex + 1 <= 1){
	return;
    }
    my $midIndex = int( ($startIndex + $endIndex) / 2);
    mergesort2arrays( $keysRef, $otherRef, $startIndex, $midIndex, $tempKeysArrayRef, $tempOtherArrayRef);
    mergesort2arrays( $keysRef, $otherRef, $midIndex + 1, $endIndex, $tempKeysArrayRef, $tempOtherArrayRef);
    merge2arrays($keysRef, $otherRef, $startIndex, $endIndex, $tempKeysArrayRef, $tempOtherArrayRef);
}

#
# Quicksort (in-place)
# [wikipedia]
#
# procedure quicksort(array, left, right)
#     if right > left
#         select a pivot index (e.g. pivotIndex := left)
#         pivotNewIndex := partition(array, left, right, pivotIndex)
#         quicksort(array, left, pivotNewIndex - 1)
#         quicksort(array, pivotNewIndex + 1, right)
#
# function partition(array, left, right, pivotIndex)
#     pivotValue := array[pivotIndex]
#     swap array[pivotIndex] and array[right] // Move pivot to end
#     storeIndex := left
#     for i  from  left to right - 1
#         if array[i] <= pivotValue 
#             swap array[i] and array[storeIndex]
#             storeIndex := storeIndex + 1
#     swap array[storeIndex] and array[right] // Move pivot to its final place
#     return storeIndex


sub quickSort2Arrays{
    my ($keysRef, $otherRef) = @_;

    if( $#{$keysRef} != $#{$otherRef}){
	die( "ERROR: $#{$keysRef} != $#{$otherRef}");
    }
    
    quickSort2ArraysRecursive($keysRef, $otherRef, 0, $#{$keysRef});

    if( $DEBUG){
	for( my $i = 1; $i <= $#{$keysRef}; $i++){
	    if( $keysRef->[$i] < $keysRef->[$i - 1]){
		die( "ERROR: quickSort2Arrays did not sort the keys array (i: $i; $keysRef->[$i] < $keysRef->[$i-1]); "); 
	    }
	}
    }
}


# NOTE: uses $rightIndex as the pivotIndex

sub quickSort2ArraysRecursive{
    no warnings 'recursion'; # turn off recursion warnings in this block only ("Deep recursion on subroutine "main::quickSort2ArraysRecursive...")
    my ($keysRef, $otherRef, $leftIndex, $rightIndex) = @_;
    
    #print STDERR "quickSort2ArraysRecursive(keysRef, otherRef, $leftIndex, $rightIndex)\n";
    
    if( $rightIndex <= $leftIndex){
	# base case
	#print STDERR "\tbase case: $rightIndex <= $leftIndex\n";
	return;
    }

    #
    # partition procedure
    #
    my $pivotValue = $keysRef->[$rightIndex];
    #print "keysRef size: $#{$keysRef}, otherRef size: $#{$otherRef}, pivotValue: $pivotValue\n";
    
    my $nextSwapIndex = $leftIndex;
    
    my $tmp;
    for( my $i = $leftIndex; $i <= $rightIndex - 1; $i++){
	if( $keysRef->[$i] < $pivotValue){
	    if( $nextSwapIndex != $i){
		#print STDERR "\tswapping indices $i($keysRef->[$i]) and $nextSwapIndex($keysRef->[$nextSwapIndex])\n";
		$tmp = $keysRef->[$nextSwapIndex];
		$keysRef->[$nextSwapIndex] = $keysRef->[$i];
		$keysRef->[$i] = $tmp;

		$tmp = $otherRef->[$nextSwapIndex];
		$otherRef->[$nextSwapIndex] = $otherRef->[$i];
		$otherRef->[$i] = $tmp;
	    }

	    $nextSwapIndex++;
	}
    }

    $tmp = $keysRef->[$nextSwapIndex];
    $keysRef->[$nextSwapIndex] = $keysRef->[$rightIndex];
    $keysRef->[$rightIndex] = $tmp;

    $tmp = $otherRef->[$nextSwapIndex];
    $otherRef->[$nextSwapIndex] = $otherRef->[$rightIndex];
    $otherRef->[$rightIndex] = $tmp;

    
    #print STDERR "\tpivot index $nextSwapIndex($keysRef->[$nextSwapIndex])\n";
    # END partition procedure

    quickSort2ArraysRecursive($keysRef, $otherRef, $leftIndex, $nextSwapIndex - 1);
    quickSort2ArraysRecursive($keysRef, $otherRef, $nextSwapIndex + 1, $rightIndex);
} # END quickSort2Arrays()

1;
