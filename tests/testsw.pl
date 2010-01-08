#!/usr/bin/perl

# perl script to test code

use Cwd;

sub Help{
  print"perl script to test mdmol and mdmol_mpi\n";
  print" -a All\n";
  print" -d Die on any error\n";
  print" -s Serial tests\n";
  print" -p Parallel tests\n";
  print" -# tests #\n";
  print" -c don't compile\n";
  print" -v Verbose\n";
  exit(1);
}

sub getostype{
  my $ostype;
  $ostype = "?";
  open(IN, "printenv |");
  while (<IN>){
    if (/OSTYPE/){
      $ostype = (split(/=/))[1];
      chomp $ostype;
      last;
    }
  }
  close(IN);
  return $ostype;
}

sub sysran {
  my $rc = $_[0];
  if($rc==0) {
    # print "successful system call \"",$_[1],"\"\n";
  } elsif ($rc == 256) {
    die "command failed: exit($rc)\n";
  } else {
    print "\"$_[1]\" ran, but exited with signal ",$rc,"\n";
    if($Error){die;}
  }
}

sub stripls {  # strip the last $_[0] lines out of $_[1]
  my $pattern = $_[0];
  my $file = $_[1];
  my $ofile = $file."tmp";
  open(OFILE,">$ofile");
  open(FILE,"<$file");
  while (<FILE>) {
    if ( /$pattern/ ) {} #print "Matched $_";
    else {print OFILE $_ ;}
  }
  close(OFILE);
  close(FILE);
  rename $ofile,$file;
}

sub striphl {  # strip the first line out of $_[0]
  my $num = $_[0];
  my $file = $_[1];
  my $ofile = $file."tmp";

  open(OFILE,">$ofile");
  open(FILE,"<$file");
  while (<FILE>) {
    if ($. > $num) {print OFILE $_;}
    # else {print $_;}
  }
  close(OFILE);
  close(FILE);
  rename $ofile,$file;
}

sub stripll {  # strip the last $_[0] lines out of $_[1]
  my $num = $_[0];
  my $file = $_[1];
  my $ofile = $file."tmp";
  my $tnum = 1;
  open(FILE,"<$file");
  while (<FILE>) {$tnum++;} # count lines
  close(FILE);
  open(FILE,"<$file");
  open(OFILE,">$ofile");
  $num = $tnum-$num;
  while (<FILE>) {
    if ($. < $num) {print OFILE $_;}
  } 
  close(OFILE);
  close(FILE);
  rename $ofile,$file;
}

sub sysdiff {
  my $command = "diff $_[0] $_[1] > $_[2]";
  my $inum = $_[3];
  my $jnum = $_[4];
  my $knum = $_[5];
  my $rc = system($command);
  if(!($rc==0)) {
    $nerrors++;
    print "\nThere are differences between \"$_[0]\" and \"$_[1]\" i=$inum, j=$jnum, nproc=$knum\n";
    if($Error){die;}
  } else {
    unlink("$_[2]");
  # if (eof tmp){print "EOF!\n"; }else{print "junk\n";} #test for empty file
  }
}

sub configure {
  my $tdir = cwd();
  chdir "..";
  system("./configure");
  chdir $tdir;
}

sub compile {
  &configure;
  my $tdir = cwd();
  my $mk = "make -f $_[0]";
  my $command = "$mk $_[1]";
  chdir "../src";
  system("make -f $_[0] clean");
  print "Making executable $_[1] using \"$command\"\n";
  my $rc = system("$command");
  if($rc==0) {
    print "Compiled successfully\n";
  } elsif ($rc == 256) {
    die "Compilation failed.... exiting\n";
  } else {
    die "Error while compiling! ... exiting\n";
  }
  chdir $tdir;
}

# this is the start of the main routines

undef($Compile);
undef($Serial);
undef($Error);
undef($MPI);
undef($t1);
undef($t1e);
undef($t2);
undef($t2e);
undef($t3);
undef($t4);

$Compile=1;

if($#ARGV==-1){&Help;} # help the user, they forgot the syntax, otherwise...

if(($#ARGV==0) && ($ARGV[0] !~ /^-/)) {
  $Dir  = shift(@ARGV);
} else {
  foreach $arg (@ARGV)  {
    if ($arg =~ /^-/) {
      if ($arg =~ /a/)  {
	$Serial=1;$MPI=1;$Compile=1;
	$t1=1;$t1e=1;$t2=1;$t2e=1;$t3=1;$t4=1;
      }
      if ($arg =~ /s/)  {
	$Serial=1;
	$t1=1;$t1e=1;$t2=1;$t2e=1;$t3=0;$t4=0;
      }
      if ($arg =~ /p/)  {
	$MPI=1;
	$t1=0;$t1e=0;$t2=0;$t2e=0;$t3=1;$t4=1;}
      if ($arg =~ /c/)  {$Compile=0; }
      if ($arg =~ /v/)  {$Verbose=1; }
      if ($arg =~ /d/)  {$Error=1; }
      if ($arg =~ /1/)  {$Serial=1;$t1=1;$t1e=1;}
      if ($arg =~ /1e/) {$Serial=1;$t1=0;$t1e=1;$t2=0;$t2e=0;}
      if ($arg =~ /2/)  {$Serial=1;$t1=0;$t1e=0;$t2=1;$t2e=0;}
      if ($arg =~ /2e/) {$Serial=1;$t1=0;$t1e=0;$t2=0;$t2e=1;}
      if ($arg =~ /3/)  {$MPI=1;$t3=1;}
      if ($arg =~ /4/)  {$MPI=1;$t4=1;}
    }
  }
  shift;
  $Dir  = shift(@ARGV);
}

local $dir = cwd();
local $ostype = getostype();
local $mpirun = "mpirun";

$nerrs = 0;
print "Testing script, directory = $dir, ostype = $ostype\n";

$| = 1;  # fflush buffer

if($Serial){
  $makefile = "Makefile";
  $code = "./mdmol";
  if($Compile){&compile($makefile,$code);}
  system("cp -f $dir/../$code $code");
  print "\nTesting serial code=$code, Makefile=$makefile, ostype = $ostype\n";
  if($t1){$nerrs += &test1();}
  if($t1e){$nerrs += &test1e();}
  if($t2){$nerrs += &test2();}
  if($t2e){$nerrs += &test2e();}
  unlink("$code");
}

if ($MPI){
#  print "we don't have parallel code yet :-(\n";
#  die;
#  break;
  if($ostype eq "Linux"){
    print "Using linux definitions\n";
    $mpirun = "/usr/bin/mpirun";
    system("export LAMRSH=ssh");
    system("lamboot -v bhost.def.cmm");
  } elsif ($ostype eq "irix6.5"){
    print "Using sgi definitions\n";
    $mpirun = "mpirun -machinefile machines.irix";
  } else {
    print "Using default definitions\n";
    $mpirun = "mpirun";
  }
  $NPROC = 1;
  $makefile = "Makefile.mpi";
  $code = "mdmol.mpi";
  if($Compile){&compile($makefile,$code);}
  system("cp -f $dir/$code $code");
  $code = "/home/moore/md/mdmol/$code";
  print "Testing mpi code=$code, Makefile=$makefile, ostype = $ostype\n";
  if($t3){$nerrs += &test3();}
  if($t4){$nerrs += &test4();}
  unlink("$code");
}

if ($nerrs != 0){
  print "There were $nerrs errors\n";
  exit(1);
} else {
  print "There were no errors detected!\n";
}
exit(0);

# testing on 1 proc with simple 4 molecules
sub test1 {
  $nerrors = 0;
  print "\nTest 1.... \n";

  my $ofile = "tmp.in";
  my $tmpf = "lj.log";
  my $natoms = "<natoms> 4 </natoms>";
  my $systemfile = "<systemfile> lj.sys </systemfile>";
  my $paramfile = "<paramfile> lj.parm </paramfile>";
  my $nsteps = "<nsteps> 1000 </nsteps>";
  my $dt = "<timestep> 0.01 </timestep>";
  my $box = "<box> 10 </box>";
  my $coords = "<coordsfile> lj.4.coord </coordsfile>";
  my $rest = "<restartfile> lj.rest </restartfile>";
  my $traj = "<trajectoryfile> lj.traj </trajectoryfile>";
  my $eng = "<energyfile> lj.eng </energyfile>";
  my $log = "<logfile> $tmpf </logfile>";
  my $ncell = "<ncell> 3 </ncell>";
  my $rcut = "<rcut> 2.5 </rcut>";
  my $skin = "<skin> 1 </skin>";
  for($j=0;$j<4;$j++){
    my $command;
    my $rc;
    print "lj.4.*.$j .... ";
    for($i=0;$i<6;$i++){
      print ".$k";
      open (OFILE,">$ofile");
      print OFILE "<simulation>";
      print OFILE "$natoms\n$nsteps\n$dt\n$box\n$systemfile\n$paramfile\n";
      print OFILE "<neighbor> $i </neighbor>\n";
      print OFILE "<periodicity> $j </periodicity>\n";
      print OFILE "$coords\n$rest\n$traj\n$eng\n$log\n";
      print OFILE "$ncell\n$rcut\n$skin\n$kmax\n$alpha\n$nfreeze\n";
      print OFILE "</simulation>";
      close (OFILE);
      $command = "$code $ofile > /dev/null 2>&1";
      $rc = system($command);
      sysran($rc,$command);
      striphl(18,"$tmpf");
      stripll(2,"$tmpf");
      sysdiff("$tmpf","lj.4._.$j.out","diff4.$i.$j.out",$i,$j,0);
    }
    print("\n");
  }
  unlink "$tmpf";
  unlink "$ofile";
  $rc = $i*$j;
  if($nerrors==0) {print "Test 1.... Passed\n";}
  else {
    print "Test 1 failed $nerrors times out of $rc sub tests\n";
    if($Error){die;}
  }
  return $nerrors;
}

# testing on 1 proc with simple 4 molecules and electrostatic
sub test1e {
  $nerrors = 0;
  print "\nTest 1e.... \n";

  my $ofile = "tmp.in";
  my $tmpf = "lj.log";
  my $systemfile = "<systemfile> lj.sys </systemfile>";
  my $paramfile = "<paramfile> lj.parm </paramfile>";
  my $natoms = "<natoms> 4 </natoms>";
  my $nsteps = "<nsteps> 1000 </nsteps>";
  my $dt = "<timestep> 0.01 </timestep>";
  my $box = "<box> 10 </box>";
  my $coords = "<coordsfile> lj.4e.coord </coordsfile>";
  my $rest = "<restartfile> lj.rest </restartfile>";
  my $traj = "<trajectoryfile> lj.traj </trajectoryfile>";
  my $eng = "<energyfile> lj.eng </energyfile>";
  my $log = "<logfile> $tmpf </logfile>";
  my $ncell = "<ncell> 3 </ncell>";
  my $skin = "<skin> 1 </skin>";
  my $rcut = "<rcut> 5 </rcut>";
  my $kmax = "<kmax> 10 </kmax>";
  my $alpha = "<alpha> 1 </alpha>";
  my $nfreeze = "<nfreeze> 0 </nfreeze>";

  for($j=0;$j<4;$j++){
    my $command;
    my $rc;
    if($j==3){$rcut = "<rcut> 3 </rcut>";}
    print "lj.4.*.$j .... ";
    for($i=0;$i<6;$i++){
      print ".$i";
      open (OFILE,">$ofile");
      print OFILE "<simulation>";
      print OFILE "$natoms\n$nsteps\n$dt\n$box\n";
      print OFILE "$systemfile\n$paramfile\n";
      print OFILE "<neighbor> $i </neighbor>\n";
      print OFILE "<periodicity> $j </periodicity>\n";
      print OFILE "$coords\n$rest\n$traj\n$eng\n$log\n";
      print OFILE "$ncell\n$rcut\n$skin\n$kmax\n$alpha\n$nfreeze\n";
      print OFILE "</simulation>";
      close (OFILE);
      $command = "$code $ofile > /dev/null 2>&1 ";
      $rc = system($command);
      sysran($rc,$command);
      striphl(18,"$tmpf");
      stripll(2,"$tmpf");
      sysdiff("$tmpf","lj.4e._$j.out","diff4e.$i.$j.out",$i,$j,0);
    }
    print("\n");
  }
  unlink "$tmpf";
  unlink "$ofile";
  $rc = $i*$j;
  if($nerrors==0) {print "Test 1e.... Passed\n";}
  else {
    print "Test 1e failed $nerrors times out of $rc sub tests\n";
    if($Error){die;}
  }
  return $nerrors;
}

# testing on 1 proc with simple 255 LJ molecules
sub test2 {
  $nerrors = 0;
  print "\nTest 2.... \n";

  my $ofile = "tmp.in";
  my $tmpf = "lj.log";
  my $systemfile = "<systemfile> lj.sys </systemfile>";
  my $paramfile = "<paramfile> lj.parm </paramfile>";
  my $natoms = "<natoms> 255 </natoms>";
  my $nsteps = "<nsteps> 100 </nsteps>";
  my $dt = "<timestep> 0.01 </timestep>";
  my $box = "<box> 7.70307 </box>";
  my $coords = "<coordsfile> lj.255.coord </coordsfile>";
  my $rest = "<restartfile> lj.rest </restartfile>";
  my $traj = "<trajectoryfile> lj.traj </trajectoryfile>";
  my $eng = "<energyfile> lj.eng </energyfile>";
  my $log = "<logfile> $tmpf </logfile>";
  my $ncell = "<ncell> 3 </ncell>";
  my $rcut = "<rcut> 2.5 </rcut>";
  my $skin = "<skin> 1 </skin>";
  my $kmax = "<kmax> 10 </kmax>";
  my $alpha = "<alpha> 1 </alpha>";
  my $nfreeze = "<nfreeze> 0 </nfreeze>";

  for($j=0;$j<4;$j++){
    my $command;
    my $rc;
    print "lj.255.*.$j .... ";
    for($i=0;$i<6;$i++){
      print ".$i";
      open (OFILE,">$ofile");
      print OFILE "<simulation>";
      print OFILE "$natoms\n$nsteps\n$dt\n$box\n";
      print OFILE "$systemfile\n$paramfile\n";
      print OFILE "<neighbor> $i </neighbor>\n";
      print OFILE "<periodicity> $j </periodicity>\n";
      print OFILE "$coords\n$rest\n$eng\n$traj\n$log\n";
      print OFILE "$ncell\n$rcut\n$skin\n$kmax\n$alpha\n$nfreeze\n";
      print OFILE "</simulation>";
      close (OFILE);
      $command = "$code $ofile > /dev/null 2>&1 ";
      $rc = system($command);
      sysran($rc,$command);
      striphl(18,"$tmpf");
      stripll(3,"$tmpf");
      sysdiff("$tmpf","lj.255._.$j.out","diff255.$i.$j.out",$i,$j,0);
    }
    print "\n"
  }
  unlink "$tmpf";
  unlink "$ofile";
  $rc = $i*$j;
  if($nerrors==0) {print "Test 2.... Passed\n";}
  else {
    print "Test 2 failed $nerrors times out of $rc sub tests\n";
    if($Error){die;}
  }
  return $nerrors;
}

# testing on 1 proc with simple 255 LJ + charged molecules
sub test2e {
  $nerrors = 0;
  print "\nTest 2e.... \n";

  my $ofile = "tmp.in";
  my $tmpf = "lj.log";
  my $systemfile = "<systemfile> lj.sys </systemfile>";
  my $paramfile = "<paramfile> lj.parm </paramfile>";
  my $natoms = "<natoms> 255 </natoms>";
  my $nsteps = "<nsteps> 400 </nsteps>";
  my $dt = "<timestep> 0.01 </timestep>";
  my $box = "<box> 7.70307 </box>";
  my $coords = "<coordsfile> lj.255e.coord </coordsfile>";
  my $rest = "<restartfile> lj.rest </restartfile>";
  my $traj = "<trajectoryfile> lj.traj </trajectoryfile>";
  my $eng = "<energyfile> lj.eng </energyfile>";
  my $log = "<logfile> $tmpf </logfile>";
  my $ncell = "<ncell> 3 </ncell>";
  my $rcut = "<rcut> 2.5 </rcut>";
  my $skin = "<skin> 1 </skin>";
  my $kmax = "<kmax> 10 </kmax>";
  my $alpha = "<alpha> 1 </alpha>";
  my $nfreeze = "<nfreeze> 0 </nfreeze>";

  for($j=0;$j<4;$j++){
    my $command;
    my $rc;
    print "lj.255e.*.$j .... ";
    for($i=0;$i<6;$i++){
      print ".$i";
      open (OFILE,">$ofile");
      print OFILE "<simulation>";
      print OFILE "$natoms\n$nsteps\n$dt\n$box\n";
      print OFILE "$systemfile\n$paramfile\n";
      print OFILE "<neighbor> $i </neighbor>\n";
      print OFILE "<periodicity> $j </periodicity>\n";
      print OFILE "$coords\n$rest\n$eng\n$traj\n$log\n";
      print OFILE "$ncell\n$rcut\n$skin\n$kmax\n$alpha\n$nfreeze\n";
      print OFILE "</simulation>";
      close (OFILE);
      $command = "$code $ofile > /dev/null 2>&1 ";
      $rc = system($command);
      sysran($rc,$command);
      striphl(18,"$tmpf");
      stripll(3,"$tmpf");
      sysdiff("$tmpf","lj.255e._.$j.out","diff255e.$i.$j.out",$i,$j,0);
    }
    print "\n"
  }
  unlink "$tmpf";
  unlink "$ofile";
  $rc = $i*$j;
  if($nerrors==0) {print "Test 2e.... Passed\n";}
  else {
    print "Test 2e failed $nerrors times out of $rc sub tests\n";
    if($Error){die;}
  }
  return $nerrors;
}



# testing on multiple proc with simple 4 LJ molecules
sub test3 {
  $nerrors = 0;
  print "\nTest 3.... (MPI) \n";

  my $ofile = "tmp.in";
  my $tmpf = "lj.log";
  my $natoms = "<natoms> 4 </natoms>";
  my $systemfile = "<systemfile> lj.sys </systemfile>";
  my $paramfile = "<paramfile> lj.parm </paramfile>";
  my $nsteps = "<nsteps> 1000 </nsteps>";
  my $dt = "<timestep> 0.01 </timestep>";
  my $box = "<box> 10 </box>";
  my $coords = "<coordsfile> lj.4.coord </coordsfile>";
  my $rest = "<restartfile> lj.rest </restartfile>";
  my $traj = "<trajectoryfile> lj.traj </trajectoryfile>";
  my $eng = "<energyfile> lj.eng </energyfile>";
  my $log = "<logfile> $tmpf </logfile>";
  my $ncell = "<ncell> 3 </ncell>";
  my $rcut = "<rcut> 2.5 </rcut>";
  my $skin = "<skin> 1 </skin>";
  for($j=0;$j<4;$j++){
    my $command;
    my $rc;
    print "lj.4.*.$j .... ";
    for($k=1;$k<=$NPROC;$k++){
      print ".$k";
      for($i=0;$i<6;$i++){
	print ".$k";
	open (OFILE,">$ofile");
	print OFILE "<simulation>";
	print OFILE "$natoms\n$nsteps\n$dt\n$box\n$systemfile\n$paramfile\n";
	print OFILE "<neighbor> $i </neighbor>\n";
	print OFILE "<periodicity> $j </periodicity>\n";
	print OFILE "$coords\n$rest\n$traj\n$eng\n$log\n";
	print OFILE "$ncell\n$rcut\n$skin\n$kmax\n$alpha\n$nfreeze\n";
	print OFILE "</simulation>";
	close (OFILE);
	$command = "$mpirun -np $k $code $ofile > /dev/null 2>&1 ";
	$rc = system($command);
	sysran($rc,$command);
	stripls("Running","$tmpf");
	stripls("Timing","$tmpf");
	stripls("cputime","$tmpf");
	striphl(18,"$tmpf");
	stripll(1,"$tmpf");
	sysdiff("$tmpf","lj.4._.$j.out","diff4.$i.$j.out",$i,$j,0);
      }
    }
    print("\n");
  }
  unlink "$tmpf";
  unlink "$ofile";
  $rc = $i*$j*($k-1);
  if($nerrors==0) {
    print "Test 3.... Passed $rc sub test\n";}
  else {
    print "Test 3 failed $nerrors times out of $rc sub tests\n";
    if($Error){die;}
  }
  return $nerrors;
}

# testing on multiple proc with 255 LJ molecules
sub test4 {
  $nerrors = 0;
  print "Test 4.... \n";
  
  $coords = "lj.255.coord";
  $ofile = "tmp.in";
  $tmpf = "tmp.out"; 
  my $ncell = 3;
  for($j=0;$j<4;$j++){
    my $command;
    my $rc;
    print "lj.255.*.$j ....";
    for($k=1;$k<=$NPROC;$k++){
      print ".$k";
      for($i=0;$i<4;$i++){
	$head = "# 255 100 0.01 7.70307 $i $j $ncell lj.rest lj.ham lj.conf\n";
	open (OFILE,">$ofile");
	print OFILE "$head";
	close (OFILE);
	$command = "cat $coords >> $ofile";
	$rc = system($command);
	$command = "$mpirun -np $k $code $ofile > $tmpf";
	$rc = system($command);
	sysran($rc,$command);
	stripls("Running","$tmpf");
	stripls("Timing","$tmpf");
	stripls("cputime","$tmpf");
	striphl(4,"$tmpf"); # strip head
	stripll(1,"$tmpf"); # strip tail
	sysdiff("$tmpf","lj.255._.$j.out","diff255mpi.$i.$j.$k.out",$i,$j,$k);
      }
    }
    print "\n";
  }
  unlink "$tmpf";
  unlink "$ofile";
  $rc = $i*$j*($k-1);
  if($nerrors==0) {print "Test 4.... Passed $rc sub test\n";}
  else {print "Test 4 failed $nerrors times out of $rc sub tests\n";}
  return $nerrors;
}
