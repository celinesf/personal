#!/usr/bin/perl
# Script to transform vcd file into ms-like file
# It cuts the genome into Chunks of "length" bp size [first SNPs to last SNP]
# It checks that length of regions >"minbp" bp.
# Outputs 
#   - "outputfile"= ms like file
#   - "outputfile-sum"= summary inrormation for each region
#   - "outputfile_EstL_template"= the template file use to write the summary statistics 
#      to create the input filr for EstL

# The possible statistics are:
# S S_1 S_2 S_s S_ss S_sf1 S_sf2 S_sl S_sh S_f S_f1 S_f2 S_o S_o1 S_o2 F(S) F(S_1) F(S_2) F(S_s) F(S_ss) F(S_sf1) F(S_sf2) F(S_sl) F(S_sh) F_st H_b H_w1 H_w2 Snn pi pi_1 pi_2 theta_W theta_W1 theta_W2 theta_H theta_H1 theta_H2 D D_1 D_2 H H_1 H_2 D_star D_star1 D_star2 D_prime D_prime1 D_prime2 r_square r_square1 r_square2 nH nH_1 nH_2 Rm Rm_1 Rm_2
################
use strict;
use warnings;
############# CHANGE INFORMATION  HERE
my $statfile="Orang_EstL_summaries";
my $regsumfile="Orang_msfile_EstL_template";
my $outputfile="Orang_Data_EstL";
my $liststat="S_1 S_2 S_sl S_sh S_f	F_st D_1 D_2 D_star1 D_star2 Snn r_square1  r_square2 Rm_1 Rm_2"; #list of statistocs you want
#####################################


if(open(STAT,$statfile) && open(SUM,$regsumfile))
{
    my $st=readline(*STAT);
    print $st;
   
}
