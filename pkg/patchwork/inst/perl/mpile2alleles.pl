#!/usr/bin/perl
# Markus Rasmussen, Sebastian DiLorenzo (mostly borrowed @ Seqanswers)
# samtools mpileup -f <reference>.fa <tumorfile>.bam > <output>.pileup
# samtools mpileup -uf <reference>.fa <tumorfile>.bam | bcftools view -vcg - > <output>.vcf

use warnings;
use strict;

@ARGV == 2 or die "Need two input arguments.";

#$pileup=$ARGV[0];
#$vcf=$ARGV[1];

open(PILEUP,$ARGV[0]) or die "Could not read from $ARGV[0] , stopping.";
open(VCF,$ARGV[1]) or die "Could not read from $ARGV[1] , stopping.";

#For testing purposes
#open(OUT,">out.txt");

my %snps;
my %iupac;
$iupac{'A'}{'R'} = 'G';
$iupac{'A'}{'W'} = 'T';
$iupac{'A'}{'M'} = 'C';
$iupac{'A'}{'Y'} = 'C|T';
$iupac{'A'}{'K'} = 'G|T';
$iupac{'A'}{'S'} = 'G|C';
$iupac{'G'}{'R'} = 'A';
$iupac{'G'}{'S'} = 'C';
$iupac{'G'}{'K'} = 'T';
$iupac{'G'}{'Y'} = 'C|T';
$iupac{'G'}{'W'} = 'A|T';
$iupac{'G'}{'M'} = 'A|C';
$iupac{'C'}{'Y'} = 'T';
$iupac{'C'}{'S'} = 'G';
$iupac{'C'}{'M'} = 'A';
$iupac{'C'}{'R'} = 'A|G';
$iupac{'C'}{'W'} = 'A|T';
$iupac{'C'}{'K'} = 'G|T';
$iupac{'T'}{'Y'} = 'C';
$iupac{'T'}{'W'} = 'A';
$iupac{'T'}{'K'} = 'G';
$iupac{'T'}{'R'} = 'A|G';
$iupac{'T'}{'S'} = 'G|C';
$iupac{'T'}{'M'} = 'A|C';

#For testing purposes
#my $i = 1;
#perl mpile2alleles.pl SRR389821.mpileup SRR389821.vcf > test.out

while (<VCF>)
	{
	chomp;
	next if /^#/;
	my ($vchr, $vpos, $cons, $consQual) = (split /\t/, $_)[0,1,4,5];

	next if $cons eq 'N';
	next if length $cons > 1;

	while (<PILEUP>) 
		{
		chomp;
		my ($chr, $pos, $ref, $depth, $baseString) = (split /\t/, $_)[0,1,2,3,4];

		#uc upper case
		$ref = uc $ref;
		next if $ref eq '*';

		#If the position and chromosome match between pileup and vcf
		# add if smaller than
		if ($vchr eq $chr && $vpos == $pos)
			{
			my $snp = ($cons =~ m/[AGCT]/) ? $cons : $iupac{$ref}{$cons};
			$snps{$chr}{$pos}{'ref'} = $ref;
			$snps{$chr}{$pos}{'snp'} = $snp;
			$snps{$chr}{$pos}{'qual'} = $consQual;

			#For testing purposes
			#print OUT "VCF: $vchr $vpos $cons $consQual \n";
			#print OUT "PILEUP: $chr $pos $ref $depth $baseString \n";
			#print OUT "snp: $snp \n";
			#print OUT "assigned: $snps{$chr}{$pos}{'snp'} \n";
			
			#if ($i == 10)
			#{
			#	close(OUT);
			#}

			#$i++;

			my $snpCount;
			if ($snps{$chr}{$pos}{'snp'} eq 'A') {
				$snpCount = $baseString =~ tr/Aa/Aa/;
			} elsif ($snps{$chr}{$pos}{'snp'} eq 'C') {
				$snpCount = $baseString =~ tr/Cc/Cc/;
			} elsif ($snps{$chr}{$pos}{'snp'} eq 'G') {
				$snpCount = $baseString =~ tr/Gg/Gg/;
			} elsif ($snps{$chr}{$pos}{'snp'} eq 'T') {
				$snpCount = $baseString =~ tr/Tt/Tt/;
			} elsif ($snps{$chr}{$pos}{'snp'} eq 'A|C') {
				$snpCount = $baseString =~ tr/AaCc/AaCc/;
			} elsif ($snps{$chr}{$pos}{'snp'} eq 'A|G') {
				$snpCount = $baseString =~ tr/AaGg/AaGg/;
			} elsif ($snps{$chr}{$pos}{'snp'} eq 'A|T') {
				$snpCount = $baseString =~ tr/AaTt/AaTt/;
			} elsif ($snps{$chr}{$pos}{'snp'} eq 'G|C') {
				$snpCount = $baseString =~ tr/GgCc/GgCc/;
			} elsif ($snps{$chr}{$pos}{'snp'} eq 'G|T') {
				$snpCount = $baseString =~ tr/GgTt/GgTt/;
			} elsif ($snps{$chr}{$pos}{'snp'} eq 'C|T') {
				$snpCount = $baseString =~ tr/CcTt/CcTt/;
			} else {
				$snpCount = 0;
			}
			$snps{$chr}{$pos}{'depth'} = $depth;
			$snps{$chr}{$pos}{'freq'} = $snpCount;
			$snps{$chr}{$pos}{'pct'} = sprintf '%.2f', ($snpCount/$depth);

			last;
			}
			elsif ($vchr eq $chr && $vpos < $pos) {
				last;
			}
		}
	}

close PILEUP;
close VCF;


foreach my $chr (sort keys %snps) {
	foreach my $pos (sort {$a <=> $b} keys %{ $snps{$chr} }) {
		print "$chr\t$pos\t${snps{$chr}{$pos}{'ref'}} > ${snps{$chr}{$pos}{'snp'}}\t${snps{$chr}{$pos}{'qual'}}\t${snps{$chr}{$pos}{'depth'}}\t${snps{$chr}{$pos}{'freq'}}\t${snps{$chr}{$pos}{'pct'}}\n";
	}
}
