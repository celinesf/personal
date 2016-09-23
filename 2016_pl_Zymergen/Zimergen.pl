#!/usr/bin/perl -w -CSDA
use strict;
use warnings;
no autovivification qw/exists/;

## 
#  The goal is to tally up the sales results from the quarter.
#  
#  The sales report input follows a CSV format with columns like the
#  following:
#  Week    Product1  Product2  Product3 ...
#     0      568.15    180.12    513.40
#     1      581.34    312.01    480.02
#    ...
#    11      545.34    134.62    502.10
# 
#  For each week, we display the sales generated from each of N products
#  represented in the second through N+1’th columns of the CSV. The 1st
#  column indicates the week number.
#
#  Each quarter consists of 12 weeks. This method generates
#  a sales report with the following aggregate data:
# 
#  * The total value associated with each week.
#  * Identify the week with the highest sales.
#  
#  Returns:
#     a record representing a sales report containing the figures of merit 
#     described above.
##

open( F, "<file.txt");
my $nline=0;
my @header;
my %total;
my $max=0;
while(<F>){
  my @line=split("\t",$_);
  if($nline==0){
      @header=@line;
  }else{
    my $week=@line[0];
    for(my $i=1;$i<scalar(@line);$i++){
      $total{week}{$week}+= $line[$i];
      $total{products}{$header[$i]}+=$line[$i];
    }
    print $week."\t".$total{week}{$week}."\n";
    if($max < $total{week}{$week}){
        $max= $total{week}{$week};
    }
  }
}
print "$max\n";

for my $p (sort keys %{$total{products}}){
  
  my $mean=  $total{products}{$p}/12;
  print "$p\t$total{products}{$p}\t$mean\n ";

}

# Part 2: Assume the following additional requirements are added to the Sales Report:

# The record must contain:
#  -- The total sales associated with each product over the quarter.
#  -- The average weekly sales associated with each product.
#  -- Products should be indexed by product name in the report.