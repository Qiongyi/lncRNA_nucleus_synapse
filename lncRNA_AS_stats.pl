#!/usr/bin/perl -w
use warnings;
use strict;
if(@ARGV != 4) {
    print STDERR "Usage: lncRNA_AS_stats.pl lncRNA_nucleus_enriched.xls suppa_output_prefix(such as merge) outf(count_stat) outf2(.xls)\n";
  	exit(0);
}
my ($inf, $suppa, $outf, $outf2)=@ARGV;

my %hash;
my $id_index;
open(IN, $inf) or die "cannot open $inf\n";
while(<IN>){
	if($_=~/^gene_id/){
		$id_index=1;
		next;
	}elsif($_=~/^transcript_id/){
		$id_index=0;
		next;
	}
	my @info=split(/\t/, $_);
	$hash{$info[$id_index]}=1;
}
close IN;

my $count_se=0;
my $as=$suppa."_SE_strict.ioe";
open(IN, $as) or die "cannot open $as\n";
while(<IN>){
	my @info=split(/\t/, $_);
	my @ele=split(/,/, $info[3]);
	foreach my $ele (@ele){
		if(exists $hash{$ele} && $hash{$ele}!~/SE\//){
			$hash{$ele}.="SE/";
			$count_se++;
			last;
		}
	}
}
close IN;

my $count_a3=0;
$as=$suppa."_A3_strict.ioe";
open(IN, $as) or die "cannot open $as\n";
while(<IN>){
	my @info=split(/\t/, $_);
	my @ele=split(/,/, $info[3]);
	foreach my $ele (@ele){
		if(exists $hash{$ele} && $hash{$ele}!~/A3SS\//){
			$hash{$ele}.="A3SS/";
			$count_a3++;
			last;
		}
	}
}
close IN;

my $count_a5=0;
$as=$suppa."_A5_strict.ioe";
open(IN, $as) or die "cannot open $as\n";
while(<IN>){
	my @info=split(/\t/, $_);
	my @ele=split(/,/, $info[3]);
	foreach my $ele (@ele){
		if(exists $hash{$ele} && $hash{$ele}!~/A5SS\//){
			$hash{$ele}.="A5SS/";
			$count_a5++;
			last;
		}
	}
}
close IN;

my $count_mx=0;
$as=$suppa."_MX_strict.ioe";
open(IN, $as) or die "cannot open $as\n";
while(<IN>){
	my @info=split(/\t/, $_);
	my @ele=split(/,/, $info[3]);
	foreach my $ele (@ele){
		if(exists $hash{$ele} && $hash{$ele}!~/MXE\//){
			$hash{$ele}.="MXE/";
			$count_mx++;
			last;
		}
	}
}
close IN;

my $count_ri=0;
$as=$suppa."_RI_strict.ioe";
open(IN, $as) or die "cannot open $as\n";
while(<IN>){
	my @info=split(/\t/, $_);
	my @ele=split(/,/, $info[3]);
	foreach my $ele (@ele){
		if(exists $hash{$ele} && $hash{$ele}!~/RI\//){
			$hash{$ele}.="RI/";
			$count_ri++;
			last;
		}
	}
}
close IN;

my $count_af=0;
$as=$suppa."_AF_strict.ioe";
open(IN, $as) or die "cannot open $as\n";
while(<IN>){
	my @info=split(/\t/, $_);
	my @ele=split(/,/, $info[3]);
	foreach my $ele (@ele){
		if(exists $hash{$ele} && $hash{$ele}!~/AFE\//){
			$hash{$ele}.="AFE/";
			$count_af++;
			last;
		}
	}
}
close IN;

my $count_al=0;
$as=$suppa."_AL_strict.ioe";
open(IN, $as) or die "cannot open $as\n";
while(<IN>){
	my @info=split(/\t/, $_);
	my @ele=split(/,/, $info[3]);
	foreach my $ele (@ele){
		if(exists $hash{$ele} && $hash{$ele}!~/ALE\//){
			$hash{$ele}.="ALE/";
			$count_al++;
			last;
		}
	}
}
close IN;

open(OUT, ">$outf") or die "cannot open $outf\n";
print OUT "Skipped exon\tSE\t$count_se\n";
print OUT "Alternative 5\' splice sites\tA5SS\t$count_a5\n";
print OUT "Alternative 3\' splice sites\tA3SS\t$count_a3\n";
print OUT "Mutually exclusive exons\tMXE\t$count_mx\n";
print OUT "Retained intron\tRI\t$count_ri\n";
print OUT "Alternative first exon\tAFE\t$count_af\n";
print OUT "Alternative last exon\tALE\t$count_al\n";

close OUT;

open(OUT, ">$outf2") or die "cannot open $outf2\n";
open(IN, $inf) or die "cannot open $inf\n";
while(<IN>){
	chomp;
	if($_=~/^gene_id/){
		print OUT "$_\tAS_type\n";
	}else{
		my @info=split(/\t/, $_);
		if($hash{$info[1]} !~ /\//){
			print OUT "$_\tNA\n";
		}else{
			$hash{$info[1]}=~s/\/$//;
			$hash{$info[1]}=~s/^1//;

			print OUT "$_\t$hash{$info[1]}\n";
		}
	}
}
close IN;
close OUT;
