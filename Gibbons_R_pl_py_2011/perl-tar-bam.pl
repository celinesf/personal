#!/usr/bin/perl
use strict;
use warnings;

my $fname="s_";
my $endf="_qseq.txt.bam";

open(MASTER,">Tarbam.sh");
print MASTER "#!/bin/sh \n";
for(my $line=1;$line<=8; $line++)
{
        for(my $dset=1;$dset<=68; $dset++)
        {
            my $fin="$fname$line\_00";
            my $nset="$dset$endf";
          
            $fin="$fin$nset";
            print "$fin\n";
            if(open(IN,"$fin"))
            {
                close IN;
                system("tar -cf $fin.tar $fin");
                system("bzip2 $fin.tar");
                system("rm $fin")
            }
        }# loop on datasets
}# line
close MASTER
