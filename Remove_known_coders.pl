#!/usr/bin/perl -w
use strict;
use warnings;
if(@ARGV != 4) {
    print STDERR "Usage: Remove_known_coders.pl gencode.vM25.annotation.gtf Nucleus.vs.Synapse.lncRNACapture.ballgown.xls others.xls protein_coding.xls\n";
    exit(0);
}

my ($gtf, $inf, $outf, $outf2)=@ARGV;

my %hash;
open(IN, $gtf) or die $!;
while(<IN>){
	next if $_=~/^#/;
	chomp;
	my @info=split(/\t/, $_);
	if($info[2] eq "transcript" && $info[8]=~/transcript_id "([^"]+)";.+transcript_type "([^"]+)";/){
		my ($tid, $type)=($1, $2);
		if($type=~/protein_coding/){
			$hash{$tid}=2;
		}else{
			$hash{$tid}=1;
		}
	}
}
close IN;



open(IN, $inf) or die $!;
open(OUT, ">$outf") or die $!;
open(OUT2, ">$outf2") or die $!;
while(<IN>){
	if($_=~/^transcriptIDs|Gene_Id|^Row\.name/i){
		print OUT $_;
		print OUT2 $_;
	}elsif($_=~/\w/){
		my @info=split(/\t/, $_);
		if(exists $hash{$info[0]} && $hash{$info[0]}==2){
			print OUT2 $_;
		}else{
			print OUT $_;
		}
	}
}
close IN;
close OUT;
close OUT2;
