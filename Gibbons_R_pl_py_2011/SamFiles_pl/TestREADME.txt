File ex1.fa contains two sequences cut from the human genome
build36. They were exatracted with command:

  samtools faidx human_b36.fa 2:2043966-2045540 20:67967-69550

Sequence names were changed manually for simplicity. File ex1.sam.gz
contains MAQ alignments exatracted with:

  (samtools view NA18507_maq.bam 2:2044001-2045500;
   samtools view NA18507_maq.bam 20:68001-69500)

and processed with `samtools fixmate' to make it self-consistent as a
standalone alignment.

To try samtools, you may run the following commands:

  samtools faidx ex1.fa                 # index the reference FASTA
  
./samtools import nomascus2.fasta.fai s_8_0021.sam-sam1  s_8_0021.sam-bam 
./samtools sort s_8_0021.sam-bam s_2_001_sort
  ./samtools index s_8_0021.sort.bam
./samtools tview s_2_001.sort.bam nomascus2.fasta





./samtools view -bS -o s_8_0021_qseq.txt.bam s_8_0021_qseq.txt.sam
./samtools sort s_8_0021_qseq.txt.bam s_8_0021.sort
./samtools pileup -cv -f nomascus2.fasta s_8_0021.sort.bam


 ./samtools view -bS s_8_0021_qseq.txt.sam | ./samtools sort - s  | ./samtools  mpileup -f nomascus2.fasta s.bam  





------------------
nohup ./samtools merge out-Johannes.bam s_1*.bam s_2*.bam >rm12 &
nohup ./samtools merge out-Vok.bam s_3*.bam s_4*.bam >rm34 &
nohup ./samtools merge -f out-Asia.bam s_5*.bam s_6*.bam >rm56 &
nohup ./samtools merge out-China.bam s_8*.bam s_7*.bam >rm78 &


nohup ./samtools sort -m 1000000000 out-Johannes.bam out-Johannes.sort >rs12 &
nohup ./samtools sort -m 1000000000 out-China.bam out-China.sort   >rs78 &



nohup ./samtools sort -m 10000000 out-Asia.bam out-Asia.sort  >rs56 &
nohup ./samtools sort -m 10000000 out-Vok.bam out-Vok.sort  >rs34 &

nohup ./samtools  mpileup -f nomascus2.fasta out-Vok.sort.bam >rmpvok &
nohup ./samtools  mpileup -f nomascus2.fasta out-Asia.sort.bam >rmpasia &



---------------------------
nohup ./samtools merge -f s_Johannes.bam s_1*.sort.bam s_2*.sort.bam >rm12 &
nohup ./samtools merge -f s_Vok.bam s_3*.sort.bam s_4*.sort.bam >rm34 &
nohup ./samtools merge -f s_Asia.bam s_5*.sort.bam s_6*.sort.bam >rm56 &
nohup ./samtools merge -f s_China.bam s_8*.sort.bam s_7*.sort.bam >rm78 &



-----------

nohup  ./samtools mpileup -uf  nomascus2.fasta out-Vok.sort.bam | ./bcftools view -bvcg - > o-Vok.bcf > rbcv &
nohup  ./samtools mpileup -uf  nomascus2.fasta out-Asia.sort.bam | ./bcftools view -bvcg - > o-Asia.bcf >rbca &
nohup  ./samtools mpileup -uf  nomascus2.fasta out-Johannes.sort.bam | ./bcftools view -bvcg - > o-Johannes.bcf &
nohup  ./samtools mpileup -uf  nomascus2.fasta out-China.sort.bam | ./bcftools view -bvcg - > o-China.bcf &


./bcftools view o-Vok.bcf | ./vcfutils.pl varFilter -D100 > o-Vok.vcf  
./bcftools view o-Asia.bcf | ./vcfutils.pl varFilter -D100 > o-Asia.vcf 
./bcftools view o-Johannes.bcf | ./vcfutils.pl varFilter -D100 > o-Johannes.vcf 
./bcftools view o-China.bcf | ./vcfutils.pl varFilter -D100 > o-China.vcf 

nohup ./samtools merge out-Khao.bam s_9.bam s_10.bam >rm910 &
nohup ./samtools sort -m 1000000000 out-Khao.bam out-Khao.sort >rs910 &
nohup  ./samtools mpileup -uf  nomascus2.fasta out-Khao.sort.bam | ./bcftools view -bvcg - > o-Khao.bcf > rbck &
./bcftools view o-Khao.bcf | ./vcfutils.pl varFilter -D100 > o-Khao.vcf  


