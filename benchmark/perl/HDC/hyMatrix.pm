package HDC::hyMatrix;

use strict;

my $separator = '	';

my $FALSE = 0;
my $TRUE = 1;

my $VERBOSE = $FALSE;
my $DEBUG = $FALSE;
#$DEBUG = $TRUE;

sub readMatrix{
    my ($matrixFileName, $matrix, $headers) = @_;

    if( !defined( $matrixFileName)){
	die( "ERROR: matrix filename in hyMatrix::readMatrix() is not defined");
    }
    open( MATRIX, $matrixFileName) or die( "ERROR: Can not open matrix file \"$matrixFileName\": $! ");
    my @lineValues;
    
    my $header = <MATRIX>;
    my @localHeaders = split( /\s+/, $header);
    my @skipColumn;
    for( my $i = 0; $i <= $#localHeaders; $i++){
	if( $localHeaders[$i] =~ m/^\s*$/){
	    splice( @localHeaders, $i, 1);
	    $i--;
	    next;
	}elsif( $localHeaders[$i] =~ m/^\s*\-\-\-/ ||
		$localHeaders[$i] =~ m/^ambiguous/ ){
	    $skipColumn[$i] = 1;
	}else{
	    $skipColumn[$i] = 0;
	}
    }

    my $row = 0;
    while( <MATRIX>){
	if( $_ =~ m/^\s*$/ ||
	    $_ =~ m/^\s*\#/ ||
	    $_ =~ m/^\s*\-\-\-/ ||
	    $_ =~ m/^ambiguous/ ){
	    next;
	}
	
	if( $_ !~ m/^(\S+)\s+(.*)/){
	    die("ERROR: Expecting row label (separated by whitespace) ");
	}
	if( $#localHeaders < $row ||
	    $localHeaders[$row] ne $1){
	    print STDERR "\@localHeaders: @localHeaders\n";
	    die( "ERROR: Row header \"$1\" != \"$localHeaders[$row]\"! (row index: $row) ");
	}
	@lineValues = split( /\s+/, $2);

	for( my $i = 0; $i <= $#lineValues; $i++){
	    if( $lineValues[$i] =~ m/^\s*$/){
		die("ERROR: Blank cell in lineValues (i = $i, lineValues[$i] = $lineValues[$i]) ");
	    }
	    if( $skipColumn[$i]){
		next;
	    }
	    $matrix->[$row][$i] = $lineValues[$i];
	    #DEBUG("DEBUG: matrix->[$row($localHeaders[$row])][$i($localHeaders[$i])] = $matrix->[$row][$i]\n");
	}
	$row++;
    }
    close( MATRIX);

    for( my $i = 0; $i <= $#localHeaders; $i++){
	if( $skipColumn[$i]){
	    splice( @localHeaders, $i, 1);
	    splice( @skipColumn, $i, 1);
	    $i--;
	    next;
	}
    }
    
    # if headers were specified, re-order cells
    if( $#{$headers} >= 0){
	my @crossRef;
	crossReferenceHeaders( $headers, \@localHeaders, \@crossRef, $FALSE);

	my $newMatrix = [[]];
	for( my $row = 0; $row <= $#{$headers}; $row++){
	    for( my $col = 0; $col <= $#{$headers}; $col++){
		$newMatrix->[$row][$col] = $matrix->[$crossRef[$row]][$crossRef[$col]];
		#DEBUG("DEBUG: newMatrix->[$row($headers->[$row])][$col($headers->[$col])] = $newMatrix->[$row][$col]\n");
	    }
	}
	# $matrix = $newMatrix; # why doesn't this work
	for( my $row = 0; $row <= $#{$headers}; $row++){
	    for( my $col = 0; $col <= $#{$headers}; $col++){
		$matrix->[$row][$col] = $newMatrix->[$row][$col];
	    }
	}
    }else{
	for( my $col = 0; $col <= $#localHeaders; $col++){
	    $headers->[$col] = $localHeaders[$col];
	}
    }


    #for( my $row = 0; $row <= 2; $row++){
    #	for( my $col = 0; $col <= 2; $col++){
    #	    DEBUG( "DEBUG: matrix->[$row($headers->[$row])][$col($headers->[$col])] = $matrix->[$row][$col]\n");
    #	}
    #}
    
    return $matrix
}


sub matrixNormalizer{
    my ($matrix) = @_;
    #
    # Get the inner product of the matrix with itself (the lower triangle that is)
    #
    my $innerProduct = 0.0;
    my $row;
    my $col;
    
    for( $row = 0; $row <= $#{$matrix}; $row++){
	for( $col = 0; $col <= $row; $col++){
	    $innerProduct += $matrix->[$row][$col]*$matrix->[$row][$col];
	}
    }

    DEBUG("innerProduct = $innerProduct\n");

    #
    # normalize all elements of the matrix
    #
    my $scaler = sqrt($innerProduct);
    DEBUG( "scaler = $scaler\n");

    for( $row = 0; $row <= $#{$matrix}; $row++){
	for( $col = 0; $col <= $#{$matrix->[0]}; $col++){
	    $matrix->[$row][$col] = $matrix->[$row][$col] / $scaler;
	}
    }
    return $matrix;
}

sub matrixInnerProduct_lower{

    my ($matrix1, $matrix1Headers, $matrix2, $matrix2Headers) = @_;

    my $matrix1Size = $#{$matrix1Headers};
    my $matrix2Size = $#{$matrix2Headers};

    if( $matrix1Size != $matrix2Size){
	die( "ERROR: hyMatrix was asked to take the inner product of two matricies with different size (matrix1Headers: @{$matrix1Headers}; matrix2Headers: @{$matrix2Headers}) ");
    }
    
    # cross-reference headers
    
    my @crossRef;
    crossReferenceHeaders( $matrix1Headers, $matrix2Headers, \@crossRef, $FALSE);
    #my $matrix1Index = 0;
    #my $matrix2Index = 0;
    #for( $matrix1Index = 0; $matrix1Index <= $matrix1Size; $matrix1Index++){
    #	for( $matrix2Index = 0; $matrix2Index <= $matrix2Size; $matrix2Index++){
    #	    if( $matrix1Headers->[$matrix1Index] eq $matrix2Headers->[$matrix2Index]){
    #		$crossRef[$matrix1Index] = $matrix2Index;
    #		last;
    #	    }
    #	}
    #	if( $matrix2Index > $matrix2Size){
    #	    print STDERR "ALERT: Could not cross reference matrix1 header: $matrix1Headers->[$matrix1Index] with any header from the other matrix: \"@{$matrix2Headers}\"\n";
    #	}
    #}
    #
    ##if( $DEBUG){
    ##	for( my $row = 0; $row <= $#crossRef; $row++){
    ##	    print "DEBUG: crossRef[$row] ($matrix1Headers->[$row]) = $crossRef[$row] ($matrix2Headers->[$crossRef[$row]])\n";
    ##	}
    ##}

    my $prodMat = [[]];
    my $total = 0;
    for( my $row = 0; $row <= $matrix1Size; $row++){
	for( my $col = 0; $col <= $row; $col++){
		
	    $prodMat->[$row][$col] = $matrix1->[$row][$col] * $matrix2->[$crossRef[$row]][$crossRef[$col]];
	    $total += $prodMat->[$row][$col];
	}
    }
    
    ### We don't really want to normalize the sum :)
    ###$total /= ($matrix1Size * $matrix1Size)/2 + $matrix1Size/2;

    return $total;
} # END matrixInnerProduct_lower()

sub elementWiseProduct{

    my ($matrix1, $matrix1Headers, $matrix2, $matrix2Headers) = @_;

    my $matrix1Size = $#{$matrix1};
    my $matrix2Size = $#{$matrix2};

    my $matrix1SizeCols = $#{$matrix1->[0]};
    my $matrix2SizeCols = $#{$matrix2->[0]};

    if( $matrix1SizeCols != $#{$matrix1Headers}){
	die( "ERROR: The number of header columns and the number of columns for one or both of the matrices does not match up (detected in elementWiseProduct())\nmatrix1Headers: @{$matrix1Headers} ($#{$matrix1Headers}), matrix1Size: $matrix1SizeCols\n ");
    }
    if( $matrix2SizeCols != $#{$matrix2Headers}){
	die( "ERROR: The number of header columns and the number of columns for one or both of the matrices does not match up (detected in elementWiseProduct())\nmatrix2Headers: @{$matrix2Headers} ($#{$matrix2Headers}), matrix2Size: $matrix2SizeCols\n ");
    }
    
    if( $matrix1Size != $matrix2Size ||
	$matrix1SizeCols != $matrix2SizeCols){
	die( "ERROR: hyMatrix was asked to perform an element-wise product of two matrices with different size (matrix1Headers: @{$matrix1Headers}; matrix2Headers: @{$matrix2Headers}) ");
    }
    
    
    # cross-reference headers
    
    my @crossRef;
    crossReferenceHeaders( $matrix1Headers, $matrix2Headers, \@crossRef, $FALSE);

    my $prodMat = [[]];
    my $total = 0;
    for( my $row = 0; $row <= $matrix1Size; $row++){
	for( my $col = 0; $col <= $matrix1SizeCols; $col++){
	    $prodMat->[$row][$col] = $matrix1->[$row][$col] * $matrix2->[$crossRef[$row]][$crossRef[$col]];
	    #$total += $prodMat->[$row][$col];
	}
    }
    
    return $prodMat;
} # END elementWiseProduct()



sub printMatrix{
    my ($printMatrix, $matrixOutFileName, $headers, $precision) = @_;
    printMatrixMaster( $printMatrix, $matrixOutFileName, $headers, $precision, $FALSE);
}

sub printMatrix_lower{
    my ($printMatrix, $matrixOutFileName, $headers, $precision) = @_;
    printMatrixMaster( $printMatrix, $matrixOutFileName, $headers, $precision, $TRUE);
}

sub printMatrixMaster{
    my ($printMatrix, $matrixOutFileName, $headers, $precision, $printLower) = @_;

    #print "printMatrixMaster(\$printMatrix, $matrixOutFileName, \$headers, $printLower)\n";
    
    my $row;
    my $col;

    my $openedFile = $FALSE;
    my $oldFH = select;
    if( defined( $matrixOutFileName) &&
	$matrixOutFileName =~ m/^STDERR/i){
	#*OUT = *STDERR;
	select STDERR;
	#$separator = ", ";
	#print OUT "Matrix:";
	#DEBUG("Using STDERR\n");
    }elsif( defined( $matrixOutFileName) &&
	    $matrixOutFileName ne "" &&
	    $matrixOutFileName !~ m/^STDOUT/i){
	open( OUT, ">$matrixOutFileName") or die("ERROR: Can not open / create matrixOutFileName \"$matrixOutFileName\": $! ");
	select OUT;
	$openedFile = $TRUE;
	#DEBUG("Using a file \"$matrixOutFileName\"\n");
    }else{
	#*OUT = *STDOUT;
	select STDOUT;
	#$separator = ", ";
	#DEBUG("Using STDOUT\n");
    }
    #DEBUG("openedFile: $openedFile\n");
    
    if( $#{$headers} >= 0){
	print "$separator";
	print join( $separator, @{$headers});
    }
    print "\n";
    
    my $colBound;
    #print "\$#{\$printMatrix}: $#{$printMatrix}\n";
    my $printFormatStr = "";
    my $printStr = "";
    if( defined($precision)){
	$printFormatStr = "\%.${precision}f$separator";
    }
    for( $row = 0; $row <= $#{$printMatrix}; $row++){
	if( $#{$headers} >= 0){
	    print "$headers->[$row]$separator";
	}
	#print STDERR "\n\n\$#{\$printMatrix->[$row]} = $#{$printMatrix->[$row]}\n\n";

	if( $printLower == $TRUE){
	    $colBound = $row;
	}else{
	    $colBound = $#{$printMatrix->[$row]}
	}
	for( $col = 0; $col <= $colBound; $col++){
	    #print STDERR "\\n\tcol = $col\n";
 	    if( defined($precision)){
		$printStr = sprintf("$printFormatStr", $printMatrix->[$row][$col]);
		print $printStr;
	    }else{
		print "$printMatrix->[$row][$col]$separator";
	    }
	}
	print "\n";
    }
    
    if(	$openedFile == $TRUE){
	#DEBUG("openedFile: $openedFile (closing OUT)\n");
	close( OUT);
    }
    select $oldFH;
    return 1;
} # END printMatrixMaster()

#sub printMatrix_lower{
#    my ($printMatrix, $matrixOutFileName, @headers) = @_;
#
#    #print STDERR "printMatrix_lower(\$printMatrix, $matrixOutFileName, \"@headers\" ($#headers entries)\n";
#    
#    my $row;
#    my $col;
#    
#    if( defined( $matrixOutFileName) &&
#	$matrixOutFileName =~ m/^STDERR/i){
#	*OUT = *STDERR;
#	$separator = ", ";
#	print OUT "Matrix:";
#    }elsif( defined( $matrixOutFileName) &&
#	$matrixOutFileName ne ""){
#	open( OUT, ">$matrixOutFileName") or die("ERROR: Can not open / create matrixOutFileName \"$matrixOutFileName\": $! ");
#    }else{
#	*OUT = *STDOUT;
#	$separator = ", ";
#    }
#
#    my $colBound;
#
#    if( $#headers >= 0){
#	print OUT "$separator";
#	print OUT join( $separator, @headers);
#    }
#    print OUT "\n";
#    
#    for( $row = 0; $row <= $#{$printMatrix}; $row++){
#	if( $#headers >= 0){
#	    print OUT "$headers[$row]$separator";
#	}
#	#print STDERR "\n\n\$#{\$printMatrix->[$row]} = $#{$printMatrix->[$row]}\n\n";
#	for( $col = 0; $col <= $row; $col++){
#	    #print STDERR "\\n\tcol = $col\n";
# 	    print OUT "$printMatrix->[$row][$col]$separator";
#	}
#	print OUT "\n";
#    }
#    
#    if( defined( $matrixOutFileName) &&
#	$matrixOutFileName !~ m/STDERR/i &&
#	$matrixOutFileName ne ""){
#	close( OUT);
#    }
#} # END printMatrix_lower() 

#sub printMatrix_lower{
#    my ($printMatrix, $matrixOutFileName, $header) = @_;
#
#    my $row;
#    my $col;
#    
#    if( defined( $matrixOutFileName) &&
#	$matrixOutFileName ne ""){
#	open( OUT, ">$matrixOutFileName") or die("ERROR: Can not open / create matrixOutFileName \"$matrixOutFileName\": $! ");
#    }else{
#	*OUT = *STDOUT;
#    }
#
#    my $colBound;
#
#    if( defined( $header) &&
#	$header ne ""){
#	print OUT "$header";
#    }
#    for( $row = 0; $row <= $#{$printMatrix}; $row++){
#	print OUT "$header->[$row]";
#	if( $printLower == $TRUE){
#	    $colBound = $row;
#	}else{
#	    $colBound = $#{$printMatrix->[0]}
#	}
#	for( $col = 0; $col <= $colBound; $col++){
#	    printf OUT "$separator%.3f", $printMatrix->[$row][$col];
#	}
#	print OUT "\n";
#    }
#    
#    if( $matrixOutFileName ne ""){
#	close( OUT);
#    }
#}


sub reducedRowEchelonForm{
    my ($matrix) = @_;

    DEBUG( "reducedRowEchelonForm( matrix)\n");
    #printMatrix($matrix, "STDERR");
    
    # for each row 
       # if (row, row) cell is 0, switch with first row with a non-zero in this column
       # reduce row so the (row, row) cell is 1
       # add a scaled version of this row to all others to zero out the other cells in this column


    my $largestRowIndex = $#{$matrix};
    my $largestColIndex = $#{$matrix->[0]};
    DEBUG( "\tlargestRowIndex = $largestRowIndex\n");
    DEBUG( "\tlargestColIndex = $largestColIndex\n");
    my $uninformativeColumns = 0;
    my $informativeCol = 0;
    
    for(my $row = 0; $row <= $largestRowIndex; $row++){
	$informativeCol = $row + $uninformativeColumns;
	if( $informativeCol > $largestColIndex){
	    return;
	}
	DEBUG( "\trow = $row, informativeCol = $informativeCol\n");
	
	if( $matrix->[$row][$informativeCol] == 0){
	    my $otherRow = $informativeCol + 1;
	    while( $otherRow <= $largestRowIndex &&
		   $matrix->[$otherRow][$informativeCol] == 0 ){
		$otherRow++;
		DEBUG( "\totherRow = $otherRow\n");
	    }
	    if( $otherRow > $largestRowIndex){
		print STDERR "NOTE: column $informativeCol is not part of the solution\n";
		$uninformativeColumns++;
		$row--;
		next;
	    }else{
		switchRows( $matrix, $row, $otherRow);
		#printMatrix($matrix, "STDERR");
	    }
	}

	if( $matrix->[$row][$informativeCol] != 0){
	    for( my $col = $informativeCol + 1; $col <= $largestColIndex; $col++){
		$matrix->[$row][$col] /= $matrix->[$row][$informativeCol];
	    }
	}
	$matrix->[$row][$informativeCol] = 1;
	#printMatrix($matrix, "STDERR");

	for( my $otherRow = 0; $otherRow <= $largestRowIndex; $otherRow++){
	    if( $otherRow == $row){
		next;
	    }
	    multipleRowByScalarAndToOtherRow( $matrix, $row, ($matrix->[$otherRow][$informativeCol] * -1.0), $otherRow, $informativeCol);
	}
	#printMatrix($matrix, "STDERR");
    }
} # END reducedRowEchelonForm()


sub switchRows{
    my ($matrix, $row, $otherRow) = @_;
    my $tmp;

    DEBUG("\tswitchRows( matrix, $row, $otherRow)\n");
    
    my $largestColIndex = $#{$matrix->[0]};
    DEBUG( "\tlargestColIndex = $largestColIndex\n");
    
    for( my $i = 0; $i <= $largestColIndex; $i++){
	$tmp = $matrix->[$row][$i];
	$matrix->[$row][$i] = $matrix->[$otherRow][$i];
	$matrix->[$otherRow][$i] = $tmp;
    }
} # END switchRows()


# ASSUMES: that the cells in the columns < the $row^th column that are 0
sub multipleRowByScalarAndToOtherRow{
    my ($matrix, $row, $scalar, $otherRow, $informativeCol) = @_;

    if( $scalar == 0){ return; }
    
    my $largestColIndex = $#{$matrix->[0]};

    DEBUG( "\tmultipleRowByScalarAndToOtherRow( matrix, $row, $scalar, $otherRow, $informativeCol)\n");
    
    for( my $col = $row; $col <= $largestColIndex; $col++){
	$matrix->[$otherRow][$col] += $matrix->[$row][$col] * $scalar;
    }
    if( $matrix->[$otherRow][$informativeCol] != 0){
	print STDERR "ALERT: multipleRowByScalarAndToOtherRow( , $row, $scalar, $otherRow) did not produce a 0.0 at ($otherRow,$informativeCol) (= $matrix->[$otherRow][$informativeCol])!\n";
    }
} # END multipleRowByScalarAndToOtherRow()



sub getMaxValue{
    my ($matrix) = @_;

    my $row;
    my $col;
    my $numRows = $#{$matrix} + 1;
    if( $numRows <= 0){
	return "";
    }
    my $numCols = $#{$matrix->[0]} + 1;
    if( $numCols <= 0){
	return "";
    }
    my $max = $matrix->[0][0];
    
    for( $row = 0; $row < $numRows; $row++){
	for( $col = 0; $col < $numCols; $col++){
	    if( $matrix->[$row][$col] > $max){
		$max = $matrix->[$row][$col];
	    }
	}
    }
    return $max;
}

sub printMatrixStats{
    my ($matrix, $matrixOutFileName, $onlyLower) = @_;

    if( ! defined( $onlyLower)){
	$onlyLower = $FALSE;
    }

    my $openedFile = $FALSE;
    if( defined( $matrixOutFileName) &&
	$matrixOutFileName =~ m/^STDERR/i){
	*OUT = *STDERR;
    }elsif( defined( $matrixOutFileName) &&
	    $matrixOutFileName ne "" &&
	    $matrixOutFileName !~ m/^STDOUT/i){
	open( OUT, ">$matrixOutFileName") or die("ERROR: Can not open / create matrixOutFileName \"$matrixOutFileName\": $! ");
	$openedFile = $TRUE;
    }else{
	*OUT = *STDOUT;
    }

    
    my $row;
    my $col;
    my $numRows = $#{$matrix} + 1;
    if( $numRows <= 0){
	return "";
    }
    my $numCols = $#{$matrix->[0]} + 1;
    if( $numCols <= 0){
	return "";
    }
    my $min = $matrix->[0][0];
    my $max = $matrix->[0][0];
    my $sum = 0.0;
    my $counter = 0;

    my $colBound;
    for( $row = 0; $row < $numRows; $row++){
	if( $onlyLower == $TRUE){
	    $colBound = $row;
	}else{
	    $colBound = $#{$matrix->[$row]}
	}
	for( $col = 0; $col <= $colBound; $col++){
	    if( $matrix->[$row][$col] < $min){
		$min = $matrix->[$row][$col];
	    }
	    if( $matrix->[$row][$col] > $max){
		$max = $matrix->[$row][$col];
	    }
	    $sum +=  $matrix->[$row][$col];
	    $counter++;
	}
    }

    print OUT "Matrix Stats:";
    if( $onlyLower){
	print OUT " (lower triangle)";
    }
    my $average = $sum / $counter;
    #print OUT "\tmin: $min\tmax: $max\tsum: $sum\tave: $average\tcount: $counter\n";
    print OUT "(min,max, sum, ave, count)\t$min\t$max\t$sum\t$average\t$counter\n";
    
    if(	$openedFile == $TRUE){
	close( OUT);
    }
} # END printMatrixStats()



# ASSUMES a sysmetrical matrix
# Description: remove columns and row that only have zeros in them (and their headers)
sub removeRowsAndColumnsWithJustZeros{

    my ($matrix, $matrixHeaders) = @_;

    my $row;
    my $col;
    my $numRows = $#{$matrix} + 1;
    if( $numRows <= 0){
	return "";
    }
    my $numCols = $#{$matrix->[0]} + 1;
    if( $numCols <= 0){
	return "";
    }
    if( $numRows != $numCols){
	die("ERROR: hyMatrix::removeRowsAndColumnsWithJustZeros() expects a sysmetrical matrix!\n".
	    "       Found $numRows rows and $numCols columns! ");
    }

    for( $row = 0; $row < $numRows; $row++){
	my $foundNonZeroValue = 0;
	for( $col = 0; $col < $numCols; $col++){
	    if( $matrix->[$row][$col] != 0){
		$foundNonZeroValue = 1;
		last;
	    }elsif( $matrix->[$col][$row] != 0){
		die("ERROR: hyMatrix::removeRowsAndColumnsWithJustZeros() expects a sysmetrical matrix!\n".
		    "       matrix[$row][$col] = $matrix->[$row][$col], but\n".
		    "       matrix[$col][$row] = $matrix->[$col][$row]! ");
	    }
	}
	if( $foundNonZeroValue == 0){
	    splice( @{$matrixHeaders}, $row, 1);
	    splice( @{$matrix}, $row, 1);
	    $numRows--;
	    $col = $row;
	    for( my $row = 0; $row < $numRows; $row++){
		splice( @{$matrix->[$row]}, $col, 1);
	    }
	    $numCols--;
	    $row--;
	}
    }
}


sub setRowsAndCols{
    my ($matrix, $matrixHeaders, $newHeaders) = @_;

    # print STDERR "addRowsAndCols(\$matrix, \$matrixHeaders, \$newHeaders):\n".
	#"\t\$matrixHeaders: \"@{$matrixHeaders}\"\n".
	#"\t\$newHeaders:    \"@{$newHeaders}\"\n";
    
    my $row;
    my $col;
    my $numRows = $#{$matrix} + 1;
    if( $numRows <= 0){
	return "";
    }
    my $numCols = $#{$matrix->[0]} + 1;
    if( $numCols <= 0){
	return "";
    }
    if( $numRows != $numCols){
	die("ERROR: hyMatrix::removeRowsAndColumnsWithJustZeros() expects a sysmetrical matrix!\n".
	    "       Found $numRows rows and $numCols columns! ");
    }

    my $newSize = $#{$newHeaders} + 1;
    my $newMatrix = [[]];

    for( $row = 0; $row < $newSize; $row++){
	for( $col = 0; $col < $newSize; $col++){
	    $newMatrix->[$row][$col] = 0;
	}
    }
    
    #
    # look for existing header--if found, copy values
    #

    my @crossRef;
    crossReferenceHeaders( $matrixHeaders, $newHeaders, \@crossRef, $TRUE);

    for( $row = 0; $row < $numRows; $row++){
	if( $crossRef[$row] < 0){
	    next;
	}
	
	# copy existing values
	for( $col = 0; $col < $numCols; $col++){
	    if( $crossRef[$col] < 0){
		next;
	    }
	    $newMatrix->[$crossRef[$row]][$crossRef[$col]] = $matrix->[$row][$col];
	    $newMatrix->[$crossRef[$col]][$crossRef[$row]] = $matrix->[$col][$row];
	}
    }
    return $newMatrix;
}


# crossReferenceHeaders() populates the crossRef array with the indices of matrix2Headers
# e.g., mat1[0] = 'A', mat2[26]='A', then crossRef[0] = 26 (and therefore mat1[0] = mat2[crossRef[0]])
sub crossReferenceHeaders{
    my ( $matrix1Headers, $matrix2Headers, $crossRef, $paritalListOK) = @_;

    
    my $matrix1Size = $#{$matrix1Headers};
    my $matrix2Size = $#{$matrix2Headers};

    
    if( $DEBUG){
	print STDERR "crossReferenceHeaders(\$matrix1Headers, \$matrix2Headers, \$crossRef, $paritalListOK)\n";
	print STDERR "matrix1:\n";
	for( my $matrix1Index = 0; $matrix1Index <= $matrix1Size; $matrix1Index++){
	    print STDERR "$matrix1Index: $matrix1Headers->[$matrix1Index]\n";
	}
	
	print STDERR "matrix2:\n";
	for( my $matrix2Index = 0; $matrix2Index <= $matrix2Size; $matrix2Index++){
	    print STDERR "$matrix2Index: $matrix2Headers->[$matrix2Index]\n";
	}
    }
	
    
    if( ! defined( $paritalListOK)){
	$paritalListOK = $FALSE;
    }elsif( $paritalListOK == $TRUE){
	for( my $matrix1Index = 0; $matrix1Index <= $matrix1Size; $matrix1Index++){
	    $crossRef->[$matrix1Index] = -1;
	}
    }

    #print STDERR "crossReferenceHeaders(\$matrix1Headers, \$matrix2Headers, \$crossRef, $paritalListOK):\n";
    #print STDERR "                      \$matrix1Headers: @{$matrix1Headers}\n";
    #print STDERR "                      \$matrix2Headers: @{$matrix2Headers}\n";
    
    my $matrix1Index = 0;
    my $matrix2Index = 0;
    for( my $matrix1Index = 0; $matrix1Index <= $matrix1Size; $matrix1Index++){
	for( my $matrix2Index = 0; $matrix2Index <= $matrix2Size; $matrix2Index++){
	    if( lc($matrix1Headers->[$matrix1Index]) eq lc($matrix2Headers->[$matrix2Index])){
		$crossRef->[$matrix1Index] = $matrix2Index;
		last;
	    }
	}
	if( $paritalListOK == $FALSE &&
	    $matrix2Index > $matrix2Size){
	    die( "ALERT: Could not cross reference matrix1 header: $matrix1Headers->[$matrix1Index] with any header from the other matrix: \"@{$matrix2Headers}\" ");
	}
    }

    if( $DEBUG){
        print STDERR "\nDEBUG: crossRef (".($#{$crossRef} + 1).") entries\n";
    	for( my $row = 0; $row <= $#{$crossRef}; $row++){
    	    print STDERR "DEBUG: crossRef[$row] ($matrix1Headers->[$row]) = $crossRef->[$row] ($matrix2Headers->[$crossRef->[$row]])\n";
    	}
    }

} # END crossReferenceHeaders()

	    
sub matrixTest1{
    my $matrix = [[]];
    my $row;
    my $col;
    my $size = 11;
    my @headers;
    
    for( $col = 0; $col <= $size; $col++){
	$headers[$col] = $col;
    }
	
    my $counter = 1;
    for( $row = 0; $row <= $size; $row++){
	for( $col = 0; $col <= $size; $col++){
	    $matrix->[$row][$col] = $counter++;
	}
    }


    matrixNormalizer( $matrix);
    printMatrix( $matrix, "/tmp/normalized.mat", @headers);
    my $sum = matrixInnerProduct_lower($matrix, \@headers, $matrix, \@headers);
    print "TEST: $sum\n";
}




sub scaleToFirstValue{
    my ($matrix, $targetValue) = @_;

    my $row;
    my $col;
    my $numRows = $#{$matrix} + 1;
    if( $numRows <= 0){
	return "";
    }
    my $numCols = $#{$matrix->[0]} + 1;
    if( $numCols <= 0){
	return "";
    }

    if( ! defined( $targetValue) ||
	$targetValue eq ""){
	$targetValue = 100;
    }

    if( $matrix->[0][0] == 0){
	die( "ERROR: Can not scale to the first value in the matrix, since it's zero! ");
    }
    my $scaler = $targetValue / $matrix->[0][0];
    
    for( $row = 0; $row < $numRows; $row++){
	for( $col = 0; $col < $numCols; $col++){
	    $matrix->[$row][$col] *= $scaler;
	}
    }
    return $matrix;
}



sub nonDiagonal{
    my ($matrix) = @_;

    my $row;
    my $col;
    my $numRows = $#{$matrix} + 1;
    if( $numRows <= 0){
	return "";
    }
    my $numCols = $#{$matrix->[0]} + 1;
    if( $numCols <= 0){
	return "";
    }

    for( $row = 0; $row < $numRows; $row++){
	$matrix->[$row][$row] = 0.0;
    }

    return $matrix;
}

sub logOdds{

    my ($matrix, $lambda, $headers) = @_;

    if( ! defined( $lambda)){
	$lambda = 1; # no effect
    }

    my $size = $#{$matrix};
    if( $size != $#{$matrix->[0]}){
	die("ERROR: Non square matrix given to logOdds() (rows: ".($size + 1).", cols: ". ( $#{$matrix->[0]} + 1).")");
    }
    
    #
    # calculate background frequencies (f_{a})
    #
    my @colCounts;
    my $total = 0;
    my $col = 0;
    my $row = 0;
    for( $col = 0; $col <= $size; $col++){
	for( $row = 0; $row <= $size; $row++){
	    $colCounts[$col] += $matrix->[$row][$col];
	}
	$total += $colCounts[$col];
    }

    my @backgroundFrequencies;
    for( $col = 0; $col <= $size; $col++){
	$backgroundFrequencies[$col] = $colCounts[$col] / $total;
	if( defined( $headers)){
	    DEBUG( "backgroundFrequencies of $headers->[$col] = $backgroundFrequencies[$col]\n");
	}
    }

    
    #
    # calculate probabilities (p_{a,b})
    #
    my $prob = [[]];
    for( $row = 0; $row <= $size; $row++){
	for( $col = 0; $col <= $size; $col++){
	    #$prob->[$row][$col] = $matrix->[$row][$col] / $colCounts[$col];
	    $prob->[$row][$col] = $matrix->[$row][$col] / $total;
	}
    }

    if( $DEBUG ){
	print( "probabilities (p_{a,b})\n");
	printMatrix( $prob);
	print("\n");
    }

    my $scalar = 1/$lambda;

    my $subMat = [[]];
    my $ratio;
    my $backgroundFrequenciesProduct;

    for( $row = 0; $row <= $size; $row++){
	for( $col = 0; $col <= $size; $col++){
	    $backgroundFrequenciesProduct = $backgroundFrequencies[$row]*$backgroundFrequencies[$col];
	    if( $backgroundFrequenciesProduct == 0){
		$subMat->[$row][$col] = 0;
	    }else{
		$ratio = $prob->[$row][$col]/$backgroundFrequenciesProduct;
		if( $ratio == 0){
		    $subMat->[$row][$col] = 0;
		}else{
		    $subMat->[$row][$col] = $scalar * log($ratio)/ log(10); # perl's log is really ln
		    if( defined( $headers)){
			DEBUG( "s($headers->[$row],$headers->[$col]) = $scalar * log($prob->[$row][$col] / ($backgroundFrequencies[$row]*$backgroundFrequencies[$col])) = $subMat->[$row][$col]\n");
		    }
		}
	    }
	}
    }

    return $subMat;
} # END logOdds()


sub makeClustalwMatrix{
    my ($matrix, $headers) = @_;

    my $size = $#{$matrix};
    if( $size != $#{$matrix->[0]}){
	die("ERROR: Non square matrix given to logOdds() (rows: ".($size + 1).", cols: ". ( $#{$matrix->[0]} + 1).")");
    }

    $headers->[ $#{$headers} + 1] = "*";
    
    my $clustalwMatrix = [[]];
    my $i; 
    my $lowest = $matrix->[0][0];
    for( $i = 0; $i <= $size; $i++){
	my $j;
	for($j = 0; $j <= $size; $j++){
	    $clustalwMatrix->[$i][$j] = $matrix->[$i][$j];
	    if( $matrix->[$i][$j] < $lowest){
		$lowest = $matrix->[$i][$j];
	    }
	}
    }

    if( $DEBUG){ print "makeClustalwMatrix(): lowest = $lowest\n"; }
    
    for( $i = 0; $i <= $size; $i++){
	$clustalwMatrix->[$i][$size + 1] = $lowest;
	$clustalwMatrix->[$size + 1][$i] = $lowest;
    }
    $clustalwMatrix->[$size + 1][$size + 1] = 1;

    return $clustalwMatrix;
}

# based on printMatrixStats()
sub sumMatrix{
    my ($matrix, $onlyLower) = @_;

    if( ! defined( $onlyLower)){
	$onlyLower = $FALSE;
    }

    my $openedFile = $FALSE;
    
    my $row;
    my $col;
    my $numRows = $#{$matrix} + 1;
    if( $numRows <= 0){
	return "";
    }
    my $numCols = $#{$matrix->[0]} + 1;
    if( $numCols <= 0){
	return "";
    }
    my $sum = 0.0;

    my $colBound;
    for( $row = 0; $row < $numRows; $row++){
	if( $onlyLower == $TRUE){
	    $colBound = $row;
	}else{
	    $colBound = $#{$matrix->[$row]}
	}
	for( $col = 0; $col <= $colBound; $col++){
	    $sum +=  $matrix->[$row][$col];
	}
    }
    
    return $sum;
} # END sumMatrix()



sub scaleMatrix{
    my ($matrix, $scaler) = @_;

    my $row;
    my $col;
    my $numRows = $#{$matrix} + 1;
    if( $numRows <= 0){
	return "";
    }
    my $numCols = $#{$matrix->[0]} + 1;
    if( $numCols <= 0){
	return "";
    }
    my $max = $matrix->[0][0];
    
    for( $row = 0; $row < $numRows; $row++){
	for( $col = 0; $col < $numCols; $col++){
	    $matrix->[$row][$col] *= $scaler;
	}
    }
    return $matrix;
} # END scaleMatrix()

sub addMatrices{

    my ($matrix1, $matrix1Headers, $matrix2, $matrix2Headers) = @_;

    my $matrix1Size = $#{$matrix1};
    my $matrix2Size = $#{$matrix2};

    my $matrix1SizeCols = $#{$matrix1->[0]};
    my $matrix2SizeCols = $#{$matrix2->[0]};

    if( $matrix1SizeCols != $#{$matrix1Headers}){
	die( "ERROR: The number of header columns and the number of columns for one or both of the matrices does not match up (detected in addMatrices())\nmatrix1Headers: @{$matrix1Headers} ($#{$matrix1Headers}), matrix1Size: $matrix1SizeCols\n ");
    }
    if( $matrix2SizeCols != $#{$matrix2Headers}){
	die( "ERROR: The number of header columns and the number of columns for one or both of the matrices does not match up (detected in addMatrices())\nmatrix2Headers: @{$matrix2Headers} ($#{$matrix2Headers}), matrix2Size: $matrix2SizeCols\n ");
    }
    
    if( $matrix1Size != $matrix2Size ||
	$matrix1SizeCols != $matrix2SizeCols){
	die( "ERROR: hyMatrix was asked to perform an element-wise product of two matrices with different size (matrix1Headers: @{$matrix1Headers}; matrix2Headers: @{$matrix2Headers}) ");
    }
    
    
    # cross-reference headers
    
    my @crossRef;
    crossReferenceHeaders( $matrix1Headers, $matrix2Headers, \@crossRef, $FALSE);

    my $sumMat = [[]];
    my $total = 0;
    for( my $row = 0; $row <= $matrix1Size; $row++){
	for( my $col = 0; $col <= $matrix1SizeCols; $col++){
	    $sumMat->[$row][$col] = $matrix1->[$row][$col] + $matrix2->[$crossRef[$row]][$crossRef[$col]];
	    #$total += $sumMat->[$row][$col];
	}
    }
    
    return $sumMat;
} # END addMatrices()



sub DEBUG{
    my ($str, $level) = @_;
    if( ! defined( $level)){
	$level = 1;
    }
    if( $DEBUG >= $level){
	print STDERR $str;
    }
}

sub VERBOSE{
    my ($str, $level) = @_;
    if( ! defined( $level)){
	$level = 1;
    }
    if( $VERBOSE >= $level){
	print STDERR $str;
    }
}


1;
