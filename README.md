# VCFtoFinalReportForBafRegress

Convert a VCF file to the Illumina FinalReport format required as input to the contamination detection tool bafRegress

Details on bafRegress can be found here: https://genome.sph.umich.edu/wiki/BAFRegress

The input VCF file should contain a single sample genotype column

The VCF genotype FORMAT column must contain the fields `NORMX`,`NORMY`,`X`,`Y`,`BAF`

The VCF INFO column must contain fields `ALLELE_A` and `ALLELE_B` with the allele suffixed with an asterix for that allele which appears as the REF allele in column 4

### Requirements:
Should run under Python 2.7 or 3

### Usage:
parseVCFToBAFRegress.py takes as input an uncompressed VCF format file as a stream from STDIN.


    cat example.vcf | python parseVcfToBAFRegress.py > example.finalreport.txt 


In the common use case where you wish to provide a compressed vcf.gz that is tabix indexed as input, as well as filter for  PASS variants from autosomes only, make use of [bcftools](http://samtools.github.io/bcftools/bcftools.html) command to uncompress and stream a filtered VCF file:

    bcftools view -f 'PASS,.' myfile.vcf.gz 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 | python parseVcfToBAFRegress.py > myfile.finalreport.txt

### Docker:

A docker image containing, VCFToFinalReportForBafRegress, bcftools, and BAFRegress can be built using the Dockerfile:

    docker build .
    
Once built, the enclosed programs can be run as below (here executed as from within the container):

    cd /root/tools
    bcftools/bin/bcftools view -f 'PASS,.' example.vcf.gz 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 | python parseVcfToBAFRegress.py > example.finalreport.txt
    python bafRegress.py estimate --freqfile ../data/GDA.MAF.txt example.finalreport.txt > results.txt

