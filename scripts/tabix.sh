#!/bin/sh

for VCFFILE in PD*.snps.vcf.gz ; do
    tabix -f -p vcf $VCFFILE 
done

