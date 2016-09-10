#!/usr/bin/perl -w -CSDA

#Identifying Unknown DNA Quicklyclick to expand
#
#Problem
#
#The GC-content of a DNA string is given by the percentage of symbols in the string that are 'C' or 'G'. For example, the GC-content of "AGCTATAG" is 37.5%. Note that the reverse complement of any DNA string has the same GC-content.
#
#DNA strings must be labeled when they are consolidated into a database. A commonly used method of string labeling is called FASTA format. In this format, the string is introduced by a line that begins with '>', followed by some labeling information. Subsequent lines contain the string itself; the first line to begin with '>' indicates the label of the next string.
#
#In Rosalind's implementation, a string in FASTA format will be labeled by the ID "Rosalind_xxxx", where "xxxx" denotes a four-digit code between 0000 and 9999.
#
#Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).
#
#Return: The ID of the string having the highest GC-content, followed by the GC-content of that string. Rosalind allows for a default error of 0.001 in all decimal answers unless otherwise stated; please see the note on absolute error below.
#
#Sample Dataset
#
#>Rosalind_6404
#CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
#TCCCACTAATAATTCTGAGG
#>Rosalind_5959
#CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT
#ATATCCATTTGTCAGCAGACACGC
#>Rosalind_0808
#CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC
#TGGGAACCTGCGGGCAGTAGGTGGAAT
#Sample Output
#
#Rosalind_0808
#60.919540

use strict;
use IO::CaptureOutput qw/capture_exec/;
use Data::Dumper;

my $DIR   = "/Users/becquetc/Documents/Helix/GC_content";
my $input = "rosalind_gc.txt";

open( my $FIN, "<", "$DIR/$input" ) or die " Could not open $DIR/$input\n";
my $cur_name;
my $nGC;
my $length;
my %highest;
while (<$FIN>) {
	s/^\s+//;
	s/\s+$//;
	### sequence name
	if (s/^>//) {
		if ( defined $cur_name and ( not keys %highest or (100* $nGC / $length) > $highest{GC} ) ) {
			$highest{name} = $cur_name;
			$highest{GC}   = 100 * $nGC / $length;
		}
		$cur_name = $_;
		$length = $nGC = 0;
	}

	#sequence
	else {
		$length += length($_);
		my @GC;
		if ( @GC = $_ =~ /G|C/gi ) {
			$nGC += scalar(@GC);
		}
	}
}

if ( defined $cur_name and ( not keys %highest or (100* $nGC / $length) > $highest{GC} ) ) {
	$highest{name} = $cur_name;
	$highest{GC}   = 100 * $nGC / $length;
}

print "$highest{name}\n".sprintf("%.6f",$highest{GC});   
