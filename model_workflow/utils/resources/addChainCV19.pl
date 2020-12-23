#!/usr/bin/perl -w
use strict;

my $long=@ARGV;
if ($long != 2){
    print "Usage: perl $0 <inputFile> <outputFile>\n";
    exit(0);
}

my ($file,$out) = @ARGV;

my @ch = qw(A B C D E F G H I J K L M N O P Q R S T U V W X Y Z 0 1 2 3 4 5 6 7 8 9);

open OUT,">$out";

my $i = 0;
my $cont = 0;
my $resnAnt = 0;
my %done;
open FILE,"$file";
while(<FILE>){

	# Skipping Ions
	next if (/Cl-/ or /Na+/ or /K+/);

	if(/^ATOM/ or /^HETATM/){

		my $resn=substr($_,22,5);
		$resn=~s/ //g;
		$cont++ if ($resn != $resnAnt);
		my $mon=substr($_,17,3);
		$mon=~s/ //g;
		
		if ($mon =~ /5$/ and !$done{$cont} and $cont > 1){
			$i++;
			$done{$cont} = 1;
		}

		my $f = substr($_,0,21);
		my $oldChain = substr($_,21,1);
		my $s = substr($_,22);
		my $chain = $ch[$i];
		if($oldChain eq ' '){
			print OUT "$f$chain$s";
		}
		else{
			print OUT "$_";
		}
		$resnAnt = $resn;
	}
	else{
		print OUT "$_";
	}
	
	if( (/^TER/ and !$done{$cont}) or (/OXT/) ){
		$cont=0;
		$i++;
	}
}
close FILE;

close OUT;
