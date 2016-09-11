#!/usr/bin/perl -w -CSDA

#Problem
#
#Poor-quality reads can be filtered out using the FASTQ Quality Filter tool from the FASTX toolkit. A command-line version of FASTX can be downloaded for Linux or MacOS from its website. An online interface for the FASTQ Quality Filter is also available here within the Galaxy web platform.
#
#Given: A quality threshold value qq, percentage of bases pp, and set of FASTQ entries.
#
#Return: Number of reads in filtered FASTQ entries
#
#Sample Dataset
#
#20 90
#@Rosalind_0049_1
#GCAGAGACCAGTAGATGTGTTTGCGGACGGTCGGGCTCCATGTGACACAG
#+
#FD@@;C<AI?4BA:=>C<G=:AE=><A??>764A8B797@A:58:527+,
#@Rosalind_0049_2
#AATGGGGGGGGGAGACAAAATACGGCTAAGGCAGGGGTCCTTGATGTCAT
#+
#1<<65:793967<4:92568-34:.>1;2752)24')*15;1,.3*3+*!
#@Rosalind_0049_3
#ACCCCATACGGCGAGCGTCAGCATCTGATATCCTCTTTCAATCCTAGCTA
#+
#B:EI>JDB5=>DA?E6B@@CA?C;=;@@C:6D:3=@49;@87;::;;?8+
#Sample Output
#
#2

use strict;
use IO::CaptureOutput qw/capture_exec/;;

my $DIR           = "/Users/becquetc/Documents/Helix";
my $input         = "rosalind_filt.txt";
my $filter_input  = "filt_input.txt";
my $filter_output = "filt_output.txt";

open( my $FIN,  "<", "$DIR/$input" )        or die " Could not open $DIR/$input\n";
open( my $FOUT, ">", "$DIR/$filter_input" ) or die " Could not open $DIR/$filter_input\n";
my $nl = 0;
my ( $p, $q );
while (<$FIN>) {
	s/^\s+//;
	s/\s+$//;
	if ( $nl == 0 ) {
		( $q, $p ) = split( " ", $_ );
	}
	else {
		print $FOUT $_ . "\n";
	}
	$nl++;
}
close($FIN);
close($FOUT);

my $command = "$DIR/bin/fastq_quality_filter -v -Q33 -q $q -p $p -i $DIR/$filter_input -o $DIR/$filter_output ";

my ( $stdout, $stderr, $success, $exit_code ) = capture_exec($command);

my @tp1 = split( "Output: ", $stdout );
@tp1 = split( " reads.", $tp1[1] );

print "$tp1[0]\n";
