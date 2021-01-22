#!/usr/bin/perl
#  CMIP utility to prepare PDB files for COORFMT=1 or 2
#  $Date: 2006/10/31 11:25:04 $
#  $id:$
#  $Revision: 1.1.1.1 $
#
if ($#ARGV < 2) {
	print "Usage: preppdb.pl res_lib pdb_in pdb_out\n";
	exit;
};
$db=`cat $ARGV[0]`;
foreach $lin (split (/\n/,$db)) {
    next if ($lin =~ /\#/);
    $res =substr($lin, 0, 4);  
    $at = substr ($lin, 5, 4);
    $at =~ s/ //g;
    $tip = substr ($lin, 10, 4);
    $car = substr ($lin, 13);
    $C{"$res$at"}=$car;
    $T{"$res$at"}=$tip;
};
open (INFL,$ARGV[1]);
open (OUT, ">$ARGV[2]");
while (<INFL>) {
	chop;
    if ((!/ATOM/) && (!/HETATM/)) {print OUT "$_\n";next};
    $pre=substr($_, 0, 54);
    $nomat=substr($_, 12, 4);
    $nomat=substr($nomat, 1, 3).substr($nomat, 0, 1) if ($nomat =~ /^[1-9]/);
    $nomat=~s/ //g;
    $nomr=substr($_, 17, 4);
    if (!$T{"$nomr$nomat"}) {$T{"$nomr$nomat"}="?????";
    print "Unknown: $pre\n"
    }
    printf OUT "%-s%8.4f  %-s\n", $pre, $C{"$nomr$nomat"}, $T{"$nomr$nomat"};
};
