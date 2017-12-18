package HDC::hyBioScripts;

# TO-DO:

use strict;
use HDC::hyMatrix;
use HDC::Common qw(:ALL);

# # for backtranslate()
# use lib "bioperl-1.5.0";
# use Bio::Seq;
# use Bio::SeqIO;


# my $FALSE = 0;
# my $TRUE = 1;
# 
# my $VERBOSE = $FALSE;
# #my $VERBOSE = $TRUE;
# my $DEBUG = $FALSE;
# #my $DEBUG = $TRUE;


# HANDLES:
# + (non-interleaved) FASTA formated data (i.e., >taxon\nseq\n)

my $MIN_DEBUG_LEVEL = 3;  # Only display debugging information for level 3 or greater

my @alphabet = ( "A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","-");

sub readDataFile{
    my ($fileName, $taxaNames, $seqs) = @_;
    my $isAligned = $TRUE;


    if( $DEBUG >= $MIN_DEBUG_LEVEL){ print STDERR "readDataFile($fileName, @{$taxaNames}, @{$seqs})\n"; }

    
    open( FILE, "$fileName") or die("\nERROR: Can not open sequence data file \"$fileName\": $!");
    my @lines = <FILE>;
    close( FILE);

    if( $#lines < 0){
	die( "\nERROR: No sequences found in \"$fileName\"" );
    }
    
    if( $lines[0] =~ m/\s*(\d+)\s+(\d+)\s*$/){
	# PHYLIP format
	die( "ERROR: PHYLIP format NYI :(; ");
    }else{
	# assuming FASTA format
	my $counter = 0;
	my $seqLength = -1;
	for( my $i = 0; $i <= $#lines; $i++){
	    if( $lines[$i] =~ m/^\s*$/){ next; }
	    
	    my $taxonName;
	    chomp $lines[$i];
	    $taxonName = $lines[$i];
	    $taxonName =~ s/\s*>\s*//;
	    $taxonName =~ s/\s*$//;
	    $taxaNames->[$counter] = $taxonName;
	    $i++;
	    if( $i > $#lines){
		die("\nERROR: Malformed input file \"$fileName\" (line index: $i, taxon: $taxonName) ");
	    }
	    chomp $lines[$i];
	    $seqs->[$counter] = $lines[$i];
	    if( $counter > 0){
		#print STDERR "length( \$seqs->[$counter]): ". (length( $seqs->[$counter]))." ($fileName)\n";
		if( length( $seqs->[$counter]) != $seqLength){
		    $isAligned = $FALSE;
		}
	    }else{
		$seqLength = length( $seqs->[$counter]);
		#print STDERR "seqLength: $seqLength ($fileName)\n";
	    }
	    $counter++;
	}
    }
    if( $DEBUG >= $MIN_DEBUG_LEVEL){ print STDERR "\t taxaNames: @{$taxaNames}\n\t seqs: @{$seqs}\n"; }

    return $isAligned;
}

sub printDataSet{
    my ($taxaNames, $seqs, $fileName, $formatType, $taxonNameSpaces) = @_;

    my $openedFile = $FALSE;
    my $oldFH = select;
    if( defined( $fileName) &&
	$fileName =~ m/^STDERR/i){
	select STDERR;
    }elsif( defined( $fileName) &&
	    $fileName ne "" &&
	    $fileName !~ m/^STDOUT/i){
	open( OUT, ">$fileName") or die("ERROR: Can not open / create fileName \"$fileName\": $! ");
	select OUT;
	$openedFile = $TRUE;
    }else{
	select STDOUT;
    }

    
    if( defined( $formatType) &&
	$formatType =~ m/PHYLIP/i){
	# PHYLIP
	my $seq = $seqs->[0];
	$seq =~ s/\n.+$//;
	print "". ($#{$seqs} + 1). " " . length( $seq) ."\n";

	if( ! defined( $taxonNameSpaces) ||
	    $taxonNameSpaces !~ m/^\d+$/){
	    $taxonNameSpaces = 15;
	}
	
	for( my $i = 0; $i <= $#{$seqs}; $i++){
	    $seqs->[$i] =~ s/\n/\n                /g;  # that's 15 spaces
	    printf "%-${taxonNameSpaces}s $seqs->[$i]\n", $taxaNames->[$i];
	}
    }elsif( defined( $formatType) &&
	    $formatType =~ m/^TAB/i){
	# Tab-delimited
	for( my $i = 0; $i <= $#{$seqs}; $i++){
	    $seqs->[$i] =~ s/\n/\n\t/g;
	    printf "%s\t%s\n", $taxaNames->[$i], join( '\t', split( //, $seqs->[$i]));
	}
    }else{
	# FASTA
	for( my $i = 0; $i <= $#{$seqs}; $i++){
	    print ">$taxaNames->[$i]\n$seqs->[$i]\n";
	}
    }

    if(	$openedFile == $TRUE){
	close( OUT);
    }
    select $oldFH;
} # END printDataSet()

    

sub stripAlignment{
    my ($seqs) = @_;

    for( my $i = 0; $i <= $#{$seqs}; $i++){
	$seqs->[$i] =~ s/-+//g;
    }
}
    
sub applyAlignment{
    my ($taxaNamesIn, $seqsIn, $taxaNamesRef, $seqsRef, $taxaNamesOut, $seqsOut) = @_;

    #print STDERR "applyAlignment(@{$taxaNamesIn}, @{$seqsIn}, @{$taxaNamesRef}, @{$seqsRef}, @{$taxaNamesOut}, @{$seqsOut})\n";

    if( $#{$taxaNamesIn} != $#{$taxaNamesRef} ||
	$#{$seqsIn} != $#{$seqsRef}){
	die( "\nERROR: ... ");
    }

    my @crossRef;
    hyMatrix::crossReferenceHeaders($taxaNamesIn, $taxaNamesRef, \@crossRef);
    
    for( my $i = 0; $i <= $#{$seqsIn}; $i++){

	if( $taxaNamesIn->[$i] ne $taxaNamesRef->[$crossRef[$i]]){
	    print STDERR "\$taxaNamesIn->[$i] ne \$taxaNamesRef->[$crossRef[$i]]\n";
	    print STDERR "$taxaNamesIn->[$i] ne $taxaNamesRef->[$crossRef[$i]]\n";
	    die( "\nERROR: Taxa names don't match ");
	}
	$taxaNamesOut->[$i] = $taxaNamesIn->[$i];
	$seqsOut->[$i] = "";
	
	my @seqArray = split(//, $seqsIn->[$i]);
	my @seqRefArray = split(//, $seqsRef->[$crossRef[$i]]);

	my $seqArrayIndex = 0;
	my $j;
	for( $j = 0; $j <= $#seqRefArray; $j++){
	    if( $seqRefArray[$j] eq "-"){
		$seqsOut->[$i] .= "-";
	    }else{
		if( $seqArrayIndex > $#seqArray){
		    print STDERR "i: $i, j: $j ($seqArrayIndex > $#seqArray)\n";
		    print STDERR "seqsOut:     :$seqsOut->[$i]:\n";
		    print STDERR "seqArray:    :" . join( '', @seqArray) . ":\n";
		    print STDERR "seqRefArray: :" . join( '', @seqRefArray) . ":\n";
		    die( "\nERROR: ... ");
		}else{
		    $seqsOut->[$i] .= $seqArray[$seqArrayIndex++];
		}
	    }
	}
	if( $seqArrayIndex - 1 != $#seqArray){
	    print STDERR "\nERROR: $seqArrayIndex - 1 != $#seqArray\n";
	    print STDERR "i: $i ($taxaNamesOut->[$i])\n";
	    print STDERR "seqsOut:     :$seqsOut->[$i]:\n";
	    print STDERR "seqArray:    :" . join( '', @seqArray) . ":\n";
	    print STDERR "seqRefArray: :" . join( '', @seqRefArray) . ":\n";
	    die( "ERROR: $seqArrayIndex - 1 != $#seqArray ");
	}

	if( $j != $#seqRefArray + 1){
	    die( "\nERROR: ... ");
	}
    }
} # END  applyAlignment()


sub referenceCasing{
    my ($taxaNamesIn, $seqsIn, $taxaNamesRef, $seqsRef, $taxaNamesOut, $seqsOut) = @_;

    #print STDERR "applyAlignment(@{$taxaNamesIn}, @{$seqsIn}, @{$taxaNamesRef}, @{$seqsRef}, @{$taxaNamesOut}, @{$seqsOut})\n";

    if( $#{$taxaNamesIn} != $#{$taxaNamesRef} ||
	$#{$seqsIn} != $#{$seqsRef}){
	die( "\nERROR: ... ");
    }

    my @crossRef;
    hyMatrix::crossReferenceHeaders($taxaNamesIn, $taxaNamesRef, \@crossRef);

    my $matches = 0;
    my $totalChars = 0;
    
    for( my $i = 0; $i <= $#{$seqsIn}; $i++){

	if( $taxaNamesIn->[$i] ne $taxaNamesRef->[$crossRef[$i]]){
	    print STDERR "\$taxaNamesIn->[$i] ne \$taxaNamesRef->[$crossRef[$i]]\n";
	    print STDERR "$taxaNamesIn->[$i] ne $taxaNamesRef->[$crossRef[$i]]\n";
	    die( "\nERROR: Taxa names don't match ");
	}
	$taxaNamesOut->[$i] = $taxaNamesIn->[$i];
	$seqsOut->[$i] = "";
	
	my @seqArray = split(//, $seqsIn->[$i]);
	my @seqRefArray = split(//, $seqsRef->[$crossRef[$i]]);

	my $j;
	$matches = 0;
	$totalChars = 0;
	for( $j = 0; $j <= $#seqArray && $j <= $#seqRefArray; $j++){
	    if( lc($seqArray[$j]) eq lc($seqRefArray[$j])){
		#print STDERR "lc($seqArray[$j]) eq lc($seqRefArray[$j])\n";
		$seqsOut->[$i] .= uc(lc($seqArray[$j]));
		$matches++;
	    }else{
		#print STDERR "lc($seqArray[$j]) ne lc($seqRefArray[$j])\n";
		$seqsOut->[$i] .= lc($seqArray[$j]);
	    }
	    $totalChars++;
	}
	for( $j = $#seqRefArray + 1; $j <= $#seqArray; $j++){
	    $seqsOut->[$i] .= lc($seqArray[$j]);
	    $totalChars++;
	}
    }
    if( $VERBOSE){ print "referenceCasing: Reference Sum of Pairs: ". ($matches / $totalChars)."\n"; }
    return ($matches / $totalChars);
} # END  referenceCasing()

sub applyCasing{
    my ($taxaNamesIn, $seqsIn, $taxaNamesRef, $seqsRef, $taxaNamesOut, $seqsOut) = @_;

    #print STDERR "applyAlignment(@{$taxaNamesIn}, @{$seqsIn}, @{$taxaNamesRef}, @{$seqsRef}, @{$taxaNamesOut}, @{$seqsOut})\n";

    if( $#{$taxaNamesIn} != $#{$taxaNamesRef} ||
	$#{$seqsIn} != $#{$seqsRef}){
	die( "ERROR: The number of sequences and/or taxon labels do not match (seqs: $#{$seqsIn} != $#{$seqsRef}; taxa: $#{$taxaNamesIn} != $#{$taxaNamesRef}) ");
    }

    my @crossRef;
    hyMatrix::crossReferenceHeaders($taxaNamesIn, $taxaNamesRef, \@crossRef);

    for( my $i = 0; $i <= $#{$seqsIn}; $i++){

	if( $taxaNamesIn->[$i] ne $taxaNamesRef->[$crossRef[$i]]){
	    print STDERR "\$taxaNamesIn->[$i] ne \$taxaNamesRef->[$crossRef[$i]]\n";
	    print STDERR "$taxaNamesIn->[$i] ne $taxaNamesRef->[$crossRef[$i]]\n";
	    die( "\nERROR: Taxa names don't match ");
	}
	$taxaNamesOut->[$i] = $taxaNamesIn->[$i];
	$seqsOut->[$i] = "";
	
	my @seqArray = split(//, $seqsIn->[$i]);
	my @seqRefArray = split(//, $seqsRef->[$crossRef[$i]]);

	my $j;
	for( $j = 0; $j <= $#seqArray && $j <= $#seqRefArray; $j++){
	    if( $seqRefArray[$j] =~ m/^\p{IsUpper}/){
		#print STDERR "$seqRefArray[$j] is UPPERCASE\n";
		$seqsOut->[$i] .= uc($seqArray[$j]);
	    }else{
		#print STDERR "$seqRefArray[$j] is lowercase\n";
		$seqsOut->[$i] .= lc($seqArray[$j]);
	    }
	}
	for( $j = $#seqRefArray + 1; $j <= $#seqArray; $j++){
	    $seqsOut->[$i] .= lc($seqArray[$j]);
	}
    }
} # END  applyCasing()


sub convert8StateSSTo3State{
    my ($seqsIn, $seqsOut) = @_;

    # H=HGI, E=EB, C=STC;	
    my %convertHash = (
	'B' => 'E',
	'C' => 'L',
	'E' => 'E',
	'G' => 'H',
	'H' => 'H',
	'I' => 'H',
	'L' => 'L',
	'S' => 'L',
	'T' => 'L',
	'-' => '-');

    for( my $i = 0; $i <= $#{$seqsIn}; $i++){
	$seqsOut->[$i] = "";
	
	my @seqArray = split(//, $seqsIn->[$i]);

	my $j;
	for( $j = 0; $j <= $#seqArray; $j++){
	    if( ! defined($convertHash{$seqArray[$j]})){
		die( "ERROR: Unexpected secondary structure \"$seqArray[$j]\" ");
	    }
	    $seqsOut->[$i] .= $convertHash{$seqArray[$j]};
	}
    }
} # END  convert8StateSSTo3State()




# sub backtranslate{
#     my ($taxaNamesIn, $seqsIn, $taxaNamesDNA, $seqsDNA, $taxaNamesOut, $seqsOut) = @_;
# 
#     #print STDERR "backtranslate(@{$taxaNamesIn}, @{$seqsIn}, @{$taxaNamesDNA}, @{$seqsDNA}, @{$taxaNamesOut}, @{$seqsOut})\n";
# 
#     if( $#{$taxaNamesIn} != $#{$taxaNamesDNA} ||
# 	$#{$seqsIn} != $#{$seqsDNA}){
# 	die( "\nERROR: ... ");
#     }
# 
#     for( my $i = 0; $i <= $#{$seqsIn}; $i++){
# 
# 	if( $taxaNamesIn->[$i] ne $taxaNamesDNA->[$i]){
# 	    die( "\nERROR: Matching up taxa names NYI ... ");
# 	}
# 	$taxaNamesOut->[$i] = $taxaNamesIn->[$i];
# 
# 
# 	# translate the DNA (for error checking)
# 	my $inSeqIOObj;
# 	my $inSeqObj = Bio::Seq->new(-seq => $seqsDNA->[$i]);
# 	my @seqsDNATranslated = split( //, $inSeqObj->translate()->seq());
# 	if( $DEBUG > = $MIN_DEBUG_LEVEL){ print "\@seqsDNATranslated: @seqsDNATranslated\n"; }
# 
# 	if( ($#seqsDNATranslated + 1) * 3 != length($seqsDNA->[$i])){
# 	    die("ERROR: ($#seqsDNATranslated + 1) * 3 != ". (length($seqsDNA->[$i]))."");
# 	}
# 	
# 	my @seqArray = split(//, $seqsIn->[$i]);
# 	my @seqDNAArray = split(//, $seqsDNA->[$i]);
# 
# 	my $seqDNAArrayIndex = 0;
# 	my $origGaps = 0;
# 	my $j;
# 	for( $j = 0; $j <= $#seqArray; $j++){
# 	    while( $seqDNAArrayIndex <= $#seqDNAArray && $seqDNAArray[$seqDNAArrayIndex] eq "-"){
# 		$seqDNAArrayIndex++;
# 	    }
# 	    if( $seqArray[$j] eq "-"){
# 		$seqsOut->[$i] .= "---";
# 		$origGaps++;
# 	    }else{
# 		if( $seqDNAArrayIndex + 2 > $#seqDNAArray){
# 		    print STDERR "i: $i, j: $j ($seqDNAArrayIndex + 2 > $#seqDNAArray)\n";
# 		    print STDERR "seqsOut:     :@{$seqsOut}:\n";
# 		    print STDERR "seqDNAArray:    :" . join( '', @seqDNAArray) . ":\n";
# 		    print STDERR "seqArray: :" . join( '', @seqArray) . ":\n";
# 		    die( "\nERROR: ... ");
# 		}
# 
# 		#print STDERR "\$seqsDNATranslated[$j-$origGaps](".$seqsDNATranslated[$j-$origGaps].") ne \$seqArray[$j]($seqArray[$j])\n";
# 		if( $seqsDNATranslated[$j - $origGaps] ne $seqArray[$j]){
# 		    print STDERR "ALERT: Mismatched residues: translated: ".$seqsDNATranslated[$j-$origGaps]." (${seqDNAArrayIndex}-".($seqDNAArrayIndex + 2)."), orig: $seqArray[$j] ($j)\n";
# 		}
# 				
# 		$seqsOut->[$i] .= $seqDNAArray[$seqDNAArrayIndex] . $seqDNAArray[$seqDNAArrayIndex + 1] . $seqDNAArray[$seqDNAArrayIndex + 2];
# 		$seqDNAArrayIndex += 3;
# 	    }
# 	}
# 	while( $seqDNAArrayIndex <= $#seqDNAArray && $seqDNAArray[$seqDNAArrayIndex] eq "-"){
# 	    $seqDNAArrayIndex++;
# 	}
# 	if( $seqDNAArrayIndex - 1 != $#seqDNAArray){
# 	    die( "\nERROR: $seqDNAArrayIndex - 1 != $#seqDNAArray ");
# 	}
# 
# 	#if( $j != $#seqDNAArray + 1){
# 	#    die( "\nERROR: ... ");
# 	#}
#     }
# } # END backtranslate()



sub stackDataSets{
    my ($taxaNamesIn1, $seqsIn1, $taxaNamesIn2, $seqsIn2, $taxaNamesOut, $seqsOut) = @_;

    #print STDERR "combineDataSets(@{$taxaNamesIn1}, @{$seqsIn1}, @{$taxaNamesIn2}, @{$seqsIn2}, @{$taxaNamesOut}, @{$seqsOut})\n";

    if( $#{$taxaNamesIn1} != $#{$taxaNamesIn2} ||
	$#{$seqsIn1} != $#{$seqsIn2}){
	die( "\nERROR: ... ");
    }

    for( my $i = 0; $i <= $#{$seqsIn1}; $i++){

	my $j = $i;
	if( $taxaNamesIn1->[$i] ne $taxaNamesIn2->[$j]){
	    die( "\nERROR: Matching up taxa names NYI ... ");
	}

	if( length( $seqsIn1->[$i]) != length( $seqsIn2->[$j])){
	    die( "\nERROR: Mismatched sizes ... ");
	}
	$taxaNamesOut->[$i] = $taxaNamesIn1->[$j];
	$seqsOut->[$i] = $seqsIn1->[$i] . "\n" . $seqsIn2->[$j];
    }
}



sub conservationOfEachColumn{
    my ($seqsIn, $conservationChars, $conservationLevels) = @_;

    #print STDERR "conservationsLevels(@{$seqsIn}, , )\n";

    my $seqs = [[]];
    my $numSeqs = $#{$seqsIn} + 1;
    for( my $i = 0; $i < $numSeqs ; $i++){
	my @seq = split( //, $seqsIn->[$i]);
	for( my $col = 0; $col <= $#seq; $col++){
	    $seqs->[$i][$col] = $seq[$col];
	}
    }

    for( my $col = 0; $col <= $#{$seqs->[0]}; $col++){
	my %charCounts;
	#my @asciiChars;
	for( my $i = 0; $i < $numSeqs; $i++){
	    $charCounts{$seqs->[$i][$col]}++;
	}

	my $mostCommonChar = "*";
	$charCounts{$mostCommonChar} = 0;
	
	for my $char ( keys %charCounts ) {
	    if( $charCounts{$char} > $charCounts{$mostCommonChar}){
		$mostCommonChar = $char;
	    }
	}
	$conservationChars->[$col] = $mostCommonChar;
	#$conservationLevels->[$col] = $charCounts{$mostCommonChar} / $#{$seqsIn};
	if( $mostCommonChar eq "-"){
	    $conservationLevels->[$col] = 0;
	}else{
	    $conservationLevels->[$col] = $charCounts{$mostCommonChar};
	    $conservationLevels->[$col] /= $numSeqs;
	}
    }
}# END  conservationOfEachColumn()


sub columnCounts{
    my ($seqsIn, $countsMatrix, $charsRef, $reportPercentages) = @_;

    #print STDERR "columnCounts(@{$seqsIn}, \$countsMatrix, @{$charsRef}, $reportPercentages)\n";

    # parse sequences into 2D array
    my $seqs = [[]];
    my $numSeqs = $#{$seqsIn} + 1;
    for( my $i = 0; $i < $numSeqs ; $i++){
	$seqsIn->[$i] =~ s/\s//g;  # mainly to remove newline at the end 
	my @seq = split( //, $seqsIn->[$i]);
	for( my $col = 0; $col <= $#seq; $col++){
	    $seqs->[$i][$col] = $seq[$col];
	}
    }

    if( ! defined( $charsRef)){
	print STDERR "REMOVE AFTER VERIFIED THAT IT WORKS (File: ". __FILE__. " Line: ". __LINE__.")\n";
	my @chars = ( "A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","-");
        for( my $alphaIndex = 0; $alphaIndex <= $#chars; $alphaIndex++){
  		$charsRef->[$alphaIndex] = $chars[$alphaIndex];
	}
    }

    if( ! defined( $reportPercentages)){
	$reportPercentages = $FALSE;
    }


    for( my $col = 0; $col <= $#{$seqs->[0]}; $col++){
	my %charCounts;

        for( my $alphaIndex = 0; $alphaIndex <= $#{$charsRef}; $alphaIndex++){
            $charCounts{lc($charsRef->[$alphaIndex])} = 0;
            $charCounts{uc($charsRef->[$alphaIndex])} = 0;
        }

	for( my $i = 0; $i < $numSeqs; $i++){
	    $charCounts{$seqs->[$i][$col]}++;
	}

	my $total = 0;
	my $divisor = 1;
	if( $reportPercentages){
	    for my $char ( keys %charCounts ) {
		$total += $charCounts{ $char};
	    }
	    $divisor = $total;
	}

        for( my $alphaIndex = 0; $alphaIndex <= $#{$charsRef}; $alphaIndex++){
		my $char = uc($charsRef->[$alphaIndex]);
		$countsMatrix->[ $alphaIndex][$col] = $charCounts{ $char} / $divisor;
		if( lc($char) ne uc($char)){
			$char = lc( $char);
			$countsMatrix->[ $alphaIndex][$col] += $charCounts{ $char} / $divisor;
		}
	}
	#for my $char ( keys %charCounts ) {
	#    if( ! exists($alpha2Index{ lc( $char)})){
	#	die( "$char");
	#    }
	#    $countsMatrix->[ $alpha2Index{ lc( $char)} ][$col] = $charCounts{ $char} / $divisor;
	#}
    }
}# END  columnCounts()


sub consensusSequence{
    my ($seqsIn, $consensusSeq) = @_;

    #print STDERR "consensusSeq(\$seqsIn, \$consensusSeq)\n";

    my $countsMatrix = [[]];
    #my @chars = ( "A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","-");
    my @chars = ( "-","A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z");
    columnCounts( $seqsIn, $countsMatrix, \@chars);

    for( my $col = 0; $col < length($seqsIn->[0]); $col++){
	my $consensusResidueI = 0;  # the 0th index is "-", so skip over it
	for( my $residueI = $consensusResidueI + 1; $residueI <= $#chars; $residueI++){
	    if( $countsMatrix->[ $residueI][$col] > $countsMatrix->[ $consensusResidueI][$col]){
		$consensusResidueI = $residueI;
	    }
	}
	if( $countsMatrix->[ $consensusResidueI][$col] == 0){
	    if( $countsMatrix->[ 0][$col] != 0){
		$consensusResidueI = 0;
	    }else{
		die( "ERROR: Did not find any residues at column index: $col (alphabet: @chars)!");
	    }
	}
	    
	push( @{$consensusSeq}, $chars[$consensusResidueI]);
    }

    #print STDERR "END consensusSeq(\$seqsIn, $consensusSeq (\@{\$consensusSeq}:".@{$consensusSeq}.")\n";
}# END  consensusSequence()



sub x_ParseSeqsInto2DArray{
    my ($seqsIn, $seqs) = @_;
    
    # parse sequences into 2D array
    my $numSeqs = $#{$seqsIn} + 1;
    for( my $i = 0; $i < $numSeqs ; $i++){
	$seqsIn->[$i] =~ s/\s//g;  # mainly to remove newline at the end 
	my @seq = split( //, $seqsIn->[$i]);
	for( my $col = 0; $col <= $#seq; $col++){
	    $seqs->[$i][$col] = $seq[$col];
	}
    }
}

sub switchQueryGaps{

    my ($seqsIn, $seqsOut) = @_;

    my $seqs = [[]];
    x_ParseSeqsInto2DArray($seqsIn, $seqs);

    my $numSeqs = $#{$seqs} + 1;
    for( my $col = 0; $col <= $#{$seqs->[0]}; $col++){
	if( $seqs->[0][$col] eq '-'){
	    my $seqI = 1;
	    SEQI: for( $seqI = 1; $seqI < $numSeqs; $seqI++){
		if ( $seqs->[$seqI][$col] ne '-'){
		    $seqs->[0][$col] = $seqs->[$seqI][$col];
		    $seqs->[$seqI][$col] = '-';
		    last SEQI;
		}
	    }
	    if( $seqs->[0][$col] eq '-'){
		die( "ERROR: Did not find any sequences withOUT a gap in column index $col!");
	    }
	}
    }

    # convert seqs from 2D to strings
    for( my $seqI = 0; $seqI < $numSeqs; $seqI++){
	$seqsOut->[$seqI] = join('', @{$seqs->[$seqI]});
    }

} # END switchQueryGaps()


sub getResidueCounts{
    my ($seqsIn, $frequenciesHash) = @_;

    for (keys %$frequenciesHash) {
        delete $frequenciesHash->{$_};
    }
    
    my $seqs = [[]];
    x_ParseSeqsInto2DArray($seqsIn, $seqs);

    my $numSeqs = $#{$seqs} + 1;
    for( my $seqI = 0; $seqI < $numSeqs; $seqI++){
	for( my $col = 0; $col <= $#{$seqs->[$seqI]}; $col++){
	    $frequenciesHash->{$seqs->[$seqI][$col]}++;
	    #print STDERR "residue: $seqs->[$seqI][$col] (post-freq: $frequenciesHash->{$seqs->[$seqI][$col]})\n";
	}
    }
} # END getResidueCounts()

sub getAlphabet{
    return @alphabet;
}

sub getLengthCounts{
    my ($seqsIn, $lengthsHash) = @_;
    
    for (keys %$lengthsHash) {
        delete $lengthsHash->{$_};
    }

    # my $seqs = [[]];
    # x_ParseSeqsInto2DArray($seqsIn, $seqs);
    # 
    # my $numSeqs = $#{$seqs} + 1;
    # for( my $seqI = 0; $seqI < $numSeqs; $seqI++){
    # 	$lengthsHash->{ $#{$seqs->[$seqI]} + 1 }++;
    # }
    
    my $numSeqs = $#{$seqsIn} + 1;
    for( my $i = 0; $i < $numSeqs ; $i++){
	$seqsIn->[$i] =~ s/\s//g;  # mainly to remove newline at the end
	$lengthsHash->{ length($seqsIn->[$i]) }++;
    }
}

sub x_MakeFrequenciesArrayFromCounts{

    my ($countsHash, $freqsRef, $total) = @_;

    foreach my $key (sort keys %$countsHash){
	for( my $i = 0; $i < $countsHash->{$key}; $i++){
	    push( @$freqsRef, $key);
	}
    }
    #print STDERR "freqsRef: @$freqsRef\n";
}


sub makeRandomDataset{

    my ($residueCountsHash, $lengthCountsHash, $taxonNamesOutRef, $seqsOutRef, $seqCount) = @_;

    # build array of residues with frequencies based on counts
    my @allResidues;
    x_MakeFrequenciesArrayFromCounts($residueCountsHash, \@allResidues);
    
    # build array of lengths with frequencies based on counts
    my @allLengths;
    x_MakeFrequenciesArrayFromCounts($lengthCountsHash, \@allLengths);

    if( ! defined( $seqCount) || $seqCount == 0){
	$seqCount = $#allLengths + 1;
    }elsif( $seqCount < 0){
	# treat negative numbers as multiplers
	$seqCount = abs($seqCount) * ($#allLengths + 1);
    }
    
    # foreach random seq
    #   make taxon name (based on index)
    #   choose a length based on the frequencies
    #   foreach position
    #     choose a residue based on the frequencies
    for( my $seqI = 0; $seqI < $seqCount; $seqI++){
	push(@$taxonNamesOutRef, "t$seqI");
	$seqsOutRef->[$seqI] = "";
	my $randLengthsIndex = int(rand($#allLengths + 1));
	my $randLength = $allLengths[$randLengthsIndex];
	for( my $colI = 0; $colI < $randLength; $colI++){
	    my $randResidueIndex = int( rand( $#allResidues + 1) );
	    my $randResidue = $allResidues[ $randResidueIndex];
	    $seqsOutRef->[$seqI] .= $randResidue;
	}
    }
} # END makeRandomDataset()


sub sort2Arrays{
    my ($keysRef, $otherRef) = @_;

    if( $#{$keysRef} < $#{$otherRef}){
	die( "$#{$keysRef} < $#{$otherRef}!; ");
    }
    
    # # naive sorting: replace the value at index with the next lowest value
    # for( my $i = 0; $i <= $#{$keysRef}; $i++){
    # 	my $lowest = $keysRef->[$i];
    # 	my $indexOfLowest = $i;
    # 	for( my $j = $i + 1; $j <= $#{$keysRef}; $j++){
    # 	    if( $keysRef->[$j] < $lowest){
    # 		$lowest = $keysRef->[$j];
    # 		$indexOfLowest = $j;
    # 	    }
    # 	}
    # 	my $tmp = $keysRef->[$i];
    # 	$keysRef->[$i] = $keysRef->[$indexOfLowest];
    # 	$keysRef->[$indexOfLowest] = $tmp;
    # 
    # 	$tmp = $otherRef->[$i];
    # 	$otherRef->[$i] = $otherRef->[$indexOfLowest];
    # 	$otherRef->[$indexOfLowest] = $tmp;
    # }

    my @tempArrayA;
    my @tempArrayB;
    $#tempArrayA = $#{$keysRef};
    $#tempArrayB = $#{$otherRef};

    mergesort2Arrays( $keysRef, $otherRef, 0, $#{$keysRef}, \@tempArrayA, \@tempArrayB);
} # END sort2Arrays()

    
sub merge2Arrays{
    my ($keysRef, $otherRef, $startIndex, $endIndex, $tempKeysArrayRef, $tempOtherArrayRef) = @_;
    #print STDERR "merge2Arrays( \$keysRef, \$otherRef, $startIndex, $endIndex, \$tempKeysArrayRef, \$tempOtherArrayRef)\n";
    
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
} # END merge2Arrays()

sub mergesort2Arrays{
    my ($keysRef, $otherRef, $startIndex, $endIndex, $tempKeysArrayRef, $tempOtherArrayRef) = @_;
    #print STDERR "mergesort2Arrays( \$keysRef, \$otherRef, $startIndex, $endIndex, \$tempKeysArrayRef, \$tempOtherArrayRef)\n";

    # if( $#{$keysRef} != $#{$otherRef}){
    # 	die( "$#{$keysRef} != $#{$otherRef}!; ");
    # }

    if( $endIndex - $startIndex + 1 <= 1){
	return;
    }
    my $midIndex = int( ($startIndex + $endIndex) / 2);
    mergesort2Arrays( $keysRef, $otherRef, $startIndex, $midIndex, $tempKeysArrayRef, $tempOtherArrayRef);
    mergesort2Arrays( $keysRef, $otherRef, $midIndex + 1, $endIndex, $tempKeysArrayRef, $tempOtherArrayRef);
    merge2Arrays($keysRef, $otherRef, $startIndex, $endIndex, $tempKeysArrayRef, $tempOtherArrayRef);
} # END mergesort2Arrays()


# For each of the taxon in $taxaToRemoveRef, look through the taxaNames (one by one) for an exact match
# Remove the entry from both taxaNames and seqsRef.
# If no match is found, keep going.
# Future Optimization: if there's more than x (e.g., 10) taxa to remove, put the labels in a hash.
sub removeTaxa{
    my ($taxaNamesRef, $seqsRef, $taxaToRemoveRef) = @_;
    
    if( $DEBUG >= $MIN_DEBUG_LEVEL){ print STDERR "removeTaxa(@{$taxaNamesRef}, @{$seqsRef}, @{$taxaToRemoveRef})\n"; }

    my %taxaToRemoveHash;
    for( my $i = 0; $i <= $#{$taxaToRemoveRef}; $i++){
	$taxaToRemoveHash{$taxaToRemoveRef->[$i]} = 1;
    }
    if( $DEBUG >= $MIN_DEBUG_LEVEL ){ print STDERR "\t".scalar(keys %taxaToRemoveHash)." taxa to remove\n"; }
    if( $DEBUG >= $MIN_DEBUG_LEVEL ){ print STDERR "\t".scalar( $taxaNamesRef)." taxa before removing any\n"; }

    my @newTaxaNames;
    my @newSeqs;
    for( my $j = 0; $j <= $#{$taxaNamesRef}; $j++){
	my $taxon = $taxaNamesRef->[$j];
	if( ! exists( $taxaToRemoveHash{$taxon})){
	    push( @newTaxaNames, $taxaNamesRef->[$j]);
	    push( @newSeqs, $seqsRef->[$j]);
	}
    }
    @{$taxaNamesRef} = @newTaxaNames;
    @{$seqsRef} = @newSeqs;
}


1;

