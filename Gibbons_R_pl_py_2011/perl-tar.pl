#!/usr/bin/perl
use strict;
use warnings;

my $fname="s_";
my $endf="_qseq.txt";

for(my $line=1;$line<=8; $line++)
{
    for(my $pair=1;$pair<=3; $pair+=2)
    {
        for(my $dset=1;$dset<=68; $dset++)
        {
            my $fin="$fname$line\_$pair\_00";
            my $nset="$dset$endf";
            if($dset<10)
            {$nset="0$dset$endf"}
            $fin="$fin$nset";

            if(open(IN,$fin))
            {
                close IN;
                print "$fin\n";
                
                system("tar -cf $fin.tar $fin");
                system("bzip2 $fin.tar");
            }            
        }# loop on datasets
    } #pair
}# line

