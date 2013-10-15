#!/usr/bin/perl
# Script to transform vcd file into ms-like file
# It cuts the genome into Chunks of "length" bp size [first SNPs to last SNP]
# It checks that length of regions >"minbp" bp.
# Outputs 
#   - "outputfile"= ms like file
#   - "outputfile-sum"= summary inrormation for each region
#   - "outputfile_EstL_template"= the template file use to write the summary statistics 
#      to create the input filr for EstL
################
use strict;
use warnings;
############# CHANGE INFORMATION  HERE
my $datafile="Orang_Data.vcf";
my $outputfile="Orang_msfile";
my $length=50000; #size of regions
my $minbp=10;
#####################################

## files headerfo‹
my $fout="$outputfile";
open(OUT,">$fout\_EstL\_template");
print OUT "r\tName\tx\tv\tw\tn_1\tn_2\tz_s\tz_e\tm_1\tm_2\tY//\n";
close OUT;
open(OUT,">$fout-summary");
print OUT "region number\tregion Name\tstarting position\tend position\tnumber of SNPs\n";
close OUT;


my $chr=1;
my $totl=0;
my $chrl=0;# lin/ snp# in chr
my $nreg=1;
my $treg=1;
my @chrp; # position on chr 0 start, 1 end, 2 length 3= # snp
$chrp[$nreg][0]=$chrp[$nreg][1]=$chrp[$nreg][3]=0;
$chrp[$nreg][2]=1;
my @haplo; #0- positions of snp 1-10 haplotypes
open(MS,">$outputfile");

if(open(IN,$datafile))
{
    while(<IN>)
    {
        my @st=$_;
        if($totl>0)
        {
            chomp(@st);
            my @tp=split(/chr/,$st[0]);
            my @tp1=split(/\t/,$tp[1]);
            my $nchr= $tp1[0];#chr number
            if($nchr ne $chr)
            {
                if($chrp[$nreg][3]>1 && $chrp[$nreg][1]-$chrp[$nreg][0]>10)# did find snps
                { write_ms($treg,$nreg,">>$outputfile$chr",">>$fout",\@haplo,\@chrp);}

                $chr=$nchr;
                $nreg=1;
                $treg++;
                $chrp[$nreg][0]= $chrp[$nreg][1]=$tp1[1];#new start up position
                $chrp[$nreg][2]=1;
                $chrp[$nreg][3]=$chrl=0;
                # open(MS,">$outputfile$chr");
                #close MS;
            }

            my $pos=$tp1[1];#position
            if($pos-$chrp[$nreg][0]>=$length)# lenght
            {
                if($chrp[$nreg][3]>1 && $chrp[$nreg][1]-$chrp[$nreg][0]>10)# did find snps
                {
                    write_ms($treg,$nreg,">>$outputfile$chr",">>$fout",\@haplo,\@chrp);
                    $nreg++;
                    $treg++;
                }
                if($nreg>1){
                    if( $pos-$chrp[$nreg-1][1]>=$length)# dist previous snp
                    { $chrp[$nreg][0]=$pos;  }
                    else#new start up position
                    { $chrp[$nreg][1]=$chrp[$nreg][0]= $chrp[$nreg-1][1]+1; }
                }
                else
                { $chrp[$nreg][0]= $chrp[$nreg][1]=0;}

                $chrp[$nreg][2]=1;
                $chrp[$nreg][3]=$chrl=0;
            }
            $chrp[$nreg][1]=$pos;#end pos
            $chrp[$nreg][2]=$pos- $chrp[$nreg][0];#lenght
            
            my $nsnp=$chrp[$nreg][3];
            $haplo[0][$nsnp]=$pos;# record ms position
            for(my $i=0;$i<10;$i++)
            {# fill up haplotype
                @tp=split(/\//,$tp1[$i+4]); 
                $haplo[$i*2+1][$nsnp]=$tp[0]; 
                $haplo[$i*2+2][$nsnp]=$tp[1];
                if(!($tp[0] eq "0" || $tp[0] eq "1" || $tp[1] eq "0" || $tp[1] eq "1"))
                {print "PROBLEM ALLELE line $totl snp $nsnp indv $i\n"}
            }
            $chrp[$nreg][3]++;# extra snp
            
            $chrl++;
        }
        $totl++;  
    }
    close IN;
    if($chrp[$nreg][3]>1 && $chrp[$nreg][1]-$chrp[$nreg][0]>10)
    {write_ms($treg,$nreg,">>$outputfile$chr",">>$fout",\@haplo,\@chrp);}
}
close MS;

sub write_ms
{
    my $ptreg=$_[0];
    my $pnreg=$_[1];
    my $fn=$_[2];
    my $fo=$_[3];
    my @phaplo=@{$_[4]};
    my @pchrp=@{$_[5]};

    if($pchrp[$pnreg][3]>0)
    {
        print MS "// 10 10 ./ms 20 1 -r .5  $pchrp[$pnreg][2] -I 2 10 10 reg$pnreg \nsegsites: $pchrp[$pnreg][3]\npositions: ";
        
        for(my $i=0;$i<=20;$i++)# sum 3-4-5-6 13-14 17-18-19-20
        {# fill up haplotypes
            if($i==0|| ($i>=3 && $i<=6)||$i==13||$i==14||$i>=17)
            {
                for(my $j=0;$j<$pchrp[$pnreg][3];$j++)
                {
                    print MS "$phaplo[$i][$j]";
                    if($i==0){print MS" ";}
                }   print MS"\n";
            }          
        }
        for(my $i=1;$i<=16;$i++)# born 1-2 7-8-9-10-11-12 15-16
        {# fill up haplotypes
            if(($i>=7 && $i<=12)||$i==1||$i==2||$i>=15)
            {
                for(my $j=0;$j<$pchrp[$pnreg][3];$j++)
                {
                    print MS "$phaplo[$i][$j]";
                    if($i==0){print MS" ";}
                } print MS"\n"; 
            }           
        }
        print MS"\n";
    }
    #close MS;

    open(OUT,"$fo\_EstL\_template");
    print OUT "$ptreg\tchr$chr-$nreg\t1\t1\t1\t10\t10\t$pchrp[$pnreg][0]\t$pchrp[$pnreg][1]\t0\t0\t1\n";
    close OUT;
    open(OUT,"$fo-summary");
    print OUT "$ptreg\tchr$chr-$nreg\t$pchrp[$pnreg][0]\t$pchrp[$pnreg][1]\t$pchrp[$pnreg][2]\t$pchrp[$pnreg][3]\n";
    close OUT;
}
