#!/usr/bin/perl -w

# Utility to prepare PDB files to run CMIP.
# Adds C and N terminals and atom charges.

use strict;

# Input Parameters.
my $long=@ARGV;
if ($long!=2){
    print "Usage: perl $0 <pdb_input_file> <pdb_output_file>\n";
    exit(0);
}
my ($pdb,$output)=@ARGV;

# CMIP Dat Directory.
my $cmipDir = "/home/dbeltran_local/model/workflow/model_workflow/utils/resources";

# Adding C and N terminals to pdb.
&afegirCN($pdb);

# Running script preppdb to add charges to each atom of the structure.
&preppdb("$cmipDir/res.lib","tmp.pdb","$output");

# Removing temporal files.
#system ("rm tmp.pdb");

####################################################################################
#                              Subrutines                                          #
####################################################################################

# AfegirCN:
# Adding C and N Terminals to Protein Chains.
sub afegirCN {
	my $in = shift;

	my $lineRes;
	my $res;
	my $nresAnt=0;
	my $found_first=0;
	my $prime=1;
    open(OUT,">tmp.pdb");
    open(APDB,"$in");
    while(<APDB>){
		next if !(/^ATOM/);
		$res = substr($_,17,3);
		my $nres = substr($_,22,4);
		$nres += 0; 		
		if(($nres != $nresAnt)){
			if(!$prime){
				print OUT $lineRes;
				$lineRes='';
				$found_first=1;
			}
		}
		my $line = $_;		
		if(!$found_first){
			$line =~s/$res /${res}N/g;
			$prime=0;
		}
		$lineRes.=$line;
		$nresAnt = $nres;
    }	
    close APDB;
	$lineRes =~ s/$res /${res}C/g;
    print OUT $lineRes;
    close OUT;
}

sub afegirCN2 {
    # Looking for begginning and final of chains (REMARK and TER).
 	my $in = shift;
 
    open(APDB,"$in");
    my $ant;
    my $seg=0;
    my $i=0;
    my @termF;
    my @termI;
    my $prim;
    while(<APDB>){
	next if (/ZN/ or /CA/ or /MG/ or /NA/);
	last if (/^END/);
	if(/^TER/){
	    $termF[$i]="$ant";
	    $seg=1;
	    $i++;
	    next;
	}
	if(/^REMARK/){
	    $prim=1;
	    next;
	};
	if($prim){
	    $termI[$i]="$_";
	    $prim=0;
	}
	if($seg){
	    if (!($prim)){
		$termI[$i]="$_";
	    }
	    $seg=0;
 	}
	$ant=$_;
   }

    # Adding C and N Terminals to each chain.
    open(AMBPDB, ">tmp.pdb");
    seek(APDB,0,0);
    while(<APDB>){
	next if (/^REMARK/ || /^TER/ || /^END/);
        my @id = split(' ',$_);
	my $line=$_;
	my $j=0;
 	while($j<$i){
	    my @id0=split(' ',$termI[$j]);
	    my @idl=split(' ',$termF[$j]);
	    if ($id[4] == $id0[4]) {
		# Cas cadena d'un sol residu
		if($id0[4]==$idl[4]){$j++; next;} 
		$line=~s/($id0[3]) /$1N/;
	    } elsif ($id[4] == $idl[4]) {
		$line=~s/($idl[3]) /$1C/;
	    };
	    $j++;
 	}
	print AMBPDB "$line";
    };
}

# Preppdb:
# Adding charges to each atom of the structure.
# Charges are readed from file1 (usually resHet.lib).
sub preppdb {
    my ($file1,$file2,$file3)=@_;
    my $db=`cat $file1`;
    my %C;
    my %T;
    foreach my $lin (split (/\n/,$db)) {
	next if ($lin =~ /\#/);
	my $res = substr($lin, 0, 4);  
	my $at = substr($lin, 5, 4);
	$at =~s/ //g;
	my $tip = substr($lin, 10, 4);
	my $car = substr($lin, 13);
	$C{"$res$at"}=$car;
	$T{"$res$at"}=$tip;
    };
    open (INFL,$file2);
    open (OUT, ">$file3");
    while (<INFL>) {
	next if ((!/ATOM/) && (!/HETATM/));
	my $pre = substr($_, 0, 54);
	my $nomat = substr($_, 12, 4);
	if ($nomat =~ /^[1-9]/ && (substr($nomat,3) eq ' ')){
	    $pre=~s/(............)(.)(..)(.)/$1$4$3$2/;
	}
	$nomat = substr($nomat, 1, 3).substr($nomat, 0, 1) if ($nomat =~ /^[1-9]/);
	$nomat =~s/ //g;
	my $nomr = substr($_, 17, 4);
	if (!$T{"$nomr$nomat"}) 
	{
	    $T{"$nomr$nomat"} = "?????";
	    print "Unknown: $pre\n";
	}
	printf OUT "%-s%8.4f  %-s\n", $pre, $C{"$nomr$nomat"}, $T{"$nomr$nomat"};
    };
}

# AfegirCN3:
# Adding C and N Terminals to Protein Chains.
# Identifying final chain through OXT atom.
sub afegirCN3 {
	my $in = shift;

	my %pdb;
	my %cter;
	my %nter;

	# Fixing first residue NTER.
	$nter{'@'}->{1} = 1;
	
    open(APDB,"$in");
    while(<APDB>){
		next if !(/^ATOM/);

		my $res = substr($_,17,3);
		my $nres = substr($_,22,4);		
		$nres += 0;
		my $at = substr($_, 12, 5);
		$at =~s/ //g;
		my $ch = substr($_, 21, 1);
		$ch =~s/ //g;
		$ch = "@" if (!$ch);

		# Using $nres as a key, because this script is only used
		# with ptraj outputs, with no problems of residue indexes.
		$pdb{$ch}->{$nres}->{'line'} .= $_;
		$pdb{$ch}->{$nres}->{'res'} = $res;
		
		if($at =~ /OXT/){
			# Residue with OXigen Terminal (OXT) marked as CTER.
			# (Ending of a chain).			
			$cter{$ch}->{$nres} = 1;
			
			# Residue following OXigen Terminal (OXT) marked as NTER.
			# (Beggining of a chain).
			$nter{$ch}->{$nres+1} = 1;
		}
    }	
    close APDB;

    open(OUT,">tmp.pdb");

	foreach my $ch (sort keys %pdb){
	foreach my $nres (sort { $a <=> $b } keys %{$pdb{$ch}}){
	
		my $line = $pdb{$ch}->{$nres}->{'line'};
		my $res = $pdb{$ch}->{$nres}->{'res'};

		if($nter{$ch}->{$nres}){
			$line =~s/$res /${res}N/g;
		}
		elsif($cter{$ch}->{$nres}){
			$line =~s/$res /${res}C/g;
		}
	    print OUT $line;		
	}
	}

    close OUT;
}
