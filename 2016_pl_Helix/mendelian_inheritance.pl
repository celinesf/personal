#!/usr/bin/perl -w -CSDA

#Introduction to Mendelian Inheritanceclick to expand
#
#Problem
#
#
#Figure 2. The probability of any outcome (leaf) in a probability tree diagram is given by the product of probabilities from the start of the tree to the outcome. For example, the probability that X is blue and Y is blue is equal to (2/5)(1/4), or 1/10.
#Probability is the mathematical study of randomly occurring phenomena. We will model such a phenomenon with a random variable, which is simply a variable that can take a number of different distinct outcomes depending on the result of an underlying random process.
#
#For example, say that we have a bag containing 3 red balls and 2 blue balls. If we let XX represent the random variable corresponding to the color of a drawn ball, then the probability of each of the two outcomes is given by Pr(X=red)=35Pr(X=red)=35 and Pr(X=blue)=25Pr(X=blue)=25.
#
#Random variables can be combined to yield new random variables. Returning to the ball example, let YY model the color of a second ball drawn from the bag (without replacing the first ball). The probability of YY being red depends on whether the first ball was red or blue. To represent all outcomes of XX and YY, we therefore use a probability tree diagram. This branching diagram represents all possible individual probabilities for XX and YY, with outcomes at the endpoints ("leaves") of the tree. The probability of any outcome is given by the product of probabilities along the path from the beginning of the tree; see Figure 2 for an illustrative example.
#
#An event is simply a collection of outcomes. Because outcomes are distinct, the probability of an event can be written as the sum of the probabilities of its constituent outcomes. For our colored ball example, let AA be the event "YY is blue." Pr(A)Pr(A) is equal to the sum of the probabilities of two different outcomes: Pr(X=blue and Y=blue)+Pr(X=red and Y=blue)Pr(X=blue and Y=blue)+Pr(X=red and Y=blue), or 310+110=25310+110=25 (see Figure 2 above).
#
#Given: Three positive integers kk, mm, and nn, representing a population containing k+m+nk+m+n organisms: kk individuals are homozygous dominant for a factor, mm are heterozygous, and nn are homozygous recessive.
#
#Return: The probability that two randomly selected mating organisms will produce an individual possessing a dominant allele (and thus displaying the dominant phenotype). Assume that any two organisms can mate.
#
#Sample Dataset
#
#2 2 2
#Sample Output
#
#0.78333


use strict;
use IO::CaptureOutput qw/capture_exec/;
use Data::Dumper;

my $DIR   = "/Users/becquetc/Documents/Helix/mendelian_inheritance";
my $input = "rosalind_iprb.txt";

open( my $FIN, "<", "$DIR/$input" ) or die " Could not open $DIR/$input\n";
my $l = readline($FIN);

## DD, Dd, dd
my %pop_properties;
( $pop_properties{DD}{n}, $pop_properties{Dd}{n}, $pop_properties{dd}{n} ) = split( " ", $l );

my $npop = 0;
for my $geno ( sort keys %pop_properties ) {
	my @alleles = split( "", $geno );
	for my $a (@alleles) {
		$pop_properties{$geno}{alleles}{$a}++;
	}
	$npop += $pop_properties{$geno}{n};
}

### prob of recessive homozygous genotype
my $prob_dd=0;
### choose first mating partner
for my $geno1 ( "dd", "Dd" ) {
	$pop_properties{$geno1}{freq} = $pop_properties{$geno1}{n} / $npop;
	### choose second mating partner
	for my $geno2 ( "dd", "Dd" ) {
		if ( $geno1 eq $geno2 ) {
			$pop_properties{$geno1}{P2}{$geno2}{freq} = ( $pop_properties{$geno2}{n} - 1 ) / ( $npop - 1 );
		}
		else {
			$pop_properties{$geno1}{P2}{$geno2}{freq} = ( $pop_properties{$geno2}{n} ) / ( $npop - 1 );
		}
		### prob each genotypes for geno1 geno2 mating
		for my $a1 ( sort keys %{ $pop_properties{$geno1}{alleles} } ) {
			for my $a2 ( sort keys %{ $pop_properties{$geno2}{alleles} } ) {
				my $f1 = $a1 . $a2;
				if ( $f1 eq "dD" ) {
					$f1 = "Dd";
				}
				$pop_properties{$geno1}{P2}{$geno2}{geno}{$f1}++;
				$pop_properties{$geno1}{P2}{$geno2}{ngeno}++;
			}
		}
		### prob of recessive homozygous genotype
			if(exists $pop_properties{$geno1}{P2}{$geno2}{geno}{dd} ){
				$prob_dd+= 
				#prob picking geno1
                $pop_properties{$geno1}{freq}*
				#prob picking geno2
				$pop_properties{$geno1}{P2}{$geno2}{freq}*
				#prob dd for this geno1 geno2 mating
				($pop_properties{$geno1}{P2}{$geno2}{geno}{dd}/$pop_properties{$geno1}{P2}{$geno2}{ngeno});
			}
	}
}
 my $prob_DDor_Dd = sprintf("%.5f",1- $prob_dd);

print "$prob_DDor_Dd\n";

