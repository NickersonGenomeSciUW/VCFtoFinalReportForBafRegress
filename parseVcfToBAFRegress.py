#!/usr/bin/env python2.7

import sys
import gzip

## input VCF file

print("[HEADER]")
print("[Data]")
col_header = ",".join(["SNP Name","Sample Name","GC Score","Allele1 - AB", "Allele2 - AB","X","Y","X Raw","Y Raw","B Allele Freq"])
print(col_header)

for line in sys.stdin:
  row = line.rstrip().split("\t")
  if line.startswith("#CHROM"):
    sampleName = row[9]

  if not line.startswith("#"):
    row = line.rstrip().split("\t")
    chrom = row[0]
    pos = row[1]
    ID = row[2]

    infoMap = {} 
    for item in row[7].split(";"):
      #print(item)
      if("=" in item):
        k,v = item.split("=")
        infoMap[k] = v

    # Checking array ALLELE_A is the REF or ALT in the VCF
    # The VCF INFO tag ALLELE_A or ALLELE_B will have asterix after allele if it is REF 
    if("*" in infoMap["ALLELE_A"]):
      is_ALLELE_A_ref = 1
    else:
      is_ALLELE_A_ref = 0

    # set default allele to "N" for final report alleles
    fr_allele1 = "N"
    fr_allele2 = "N"

    formatTags = row[8].split(":")
    genoData = row[9].split(":")
    genoMap = dict(zip(formatTags,genoData))

    geno = genoMap["GT"]
    if(is_ALLELE_A_ref==1):
      if(geno=="0/0"):
        fr_allele1 = "A"
        fr_allele2 = "A"
      if(geno=="0/1"):
        fr_allele1 = "A"
        fr_allele2 = "B"
      if(geno=="1/1"):
        fr_allele1 = "B"
        fr_allele2 = "B"
      if(geno=="./."):
        fr_allele1 = "-"
        fr_allele2 = "-"
    else:
      if(geno=="0/0"):
        fr_allele1 = "B"
        fr_allele2 = "B"      
      if(geno=="1/0"):
        fr_allele1 = "A"
        fr_allele2 = "B"
      if(geno=="1/1"):
        fr_allele1 = "A"
        fr_allele2 = "A"
      if(geno=="./."):
        fr_allele1 = "-"
        fr_allele2 = "-"

    outrow = [ID,sampleName,genoMap["IGC"],fr_allele1,fr_allele2,genoMap["NORMX"],genoMap["NORMY"],genoMap["X"],genoMap["Y"],genoMap["BAF"]]
    print(",".join(outrow))
