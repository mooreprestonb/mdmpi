#!/usr/bin/perl

while(<>){
  if($. == 1){
    # traj keyword
  } else {
    @Fld = split(/\s+/, $_ , 11);
    if( $. == 3){
      print "# $Fld[1] $Fld[4] $Fld[7]\n";
    }
    if($Fld[0] eq "<xyz"){
      print "$Fld[2] $Fld[3] $Fld[4]\n";
    }
    if($Fld[0] eq "<hmatrix>"){
      print "$Fld[1] $Fld[2] $Fld[3] $Fld[4] $Fld[5] $Fld[6] $Fld[7] $Fld[8] $Fld[9]\n";
    }
  }
}
