#!/usr/bin/perl

$tatoms=0;
$| = 1;  # fflush stdout
while (<>) {
    if ($. == 1 ){
	@Fld = split(/\s+/, $_ , 11); # split line into at most four parts
	die if ( $Fld[0] ne '#');
	$tatoms = $Fld[1];
	$ncell = $Fld[4];
	# print "natoms = $tatoms, ncell = $ncell\n";
	$nnow = $tatoms*27;
	$ncell *= 3;
	print "# $nnow $Fld[2] $Fld[3] $ncell $Fld[5] $Fld[6] $Fld[7] $Fld[8] $Fld[9] $Fld[10]";
	$ncell /= 3;
    } else {
	$atoms[++$#atoms] = $_;
    }
}
for($l = -1;$l <= 1;$l++){
    for($k = -1;$k <= 1;$k++){
	for($j = -1;$j <= 1;$j++){
	    for($i = 0;$i < $tatoms;$i++){
		@Fld = split(/\s+/, $atoms[$i] , 4); # split line 
		$x = $Fld[0]+$j*$ncell;
		$y = $Fld[1]+$k*$ncell;
		$z = $Fld[2]+$l*$ncell;
		print "$x $y $z $Fld[3]";
	    }
	}
    }
}
