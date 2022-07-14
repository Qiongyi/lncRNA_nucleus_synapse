#!/usr/bin/perl -w
use strict;
use warnings;
if(@ARGV != 4) {
    print STDERR "Usage: Extract_known_lnRNAs.pl Nucleus.vs.Synapse.lncRNACapture.ballgown.xls nucleus_synapse_FPKM.xls nucleus_synapse_FPKM.qval.known.xls nucleus_synapse_FPKM.qval.xls\n";
    exit(0);
}

my ($bg, $inf, $outf, $outf2)=@ARGV;

my %hash;
open(IN, $bg) or die $!;
while(<IN>){
	chomp;
	my @info=split(/\t/, $_);
	$hash{$info[0]}="$info[3]\t$info[4]\t$info[5]";
}
close IN;


open(IN, $inf) or die $!;
open(OUT, ">$outf") or die $!;
open(OUT2, ">$outf2") or die $!;
while(<IN>){
	chomp;
	if($_=~/^Gene_Id|^Row\.name/i){
		print OUT "$_\tfold-change\tP-value\tQ-value\n";
		print OUT2 "$_\tfold-change\tP-value\tQ-value\n";
	}elsif($_=~/\w/){
		my @info=split(/\t/, $_);
		if(exists $hash{$info[1]}){
			print OUT2 "$_\t$hash{$info[1]}\n";
			if($info[2] ne "NA"){
				print OUT "$_\t$hash{$info[1]}\n";
			}
		}
	}
}
close IN;
close OUT;
close OUT2;
