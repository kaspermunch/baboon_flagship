import os, sys
from gwf import Workflow
from gwf.workflow import collect
from templates_vcfmerge import *

gwf = Workflow(defaults={'account': 'primatediversity'})

# paths to VCF files:
vcf_files = gwf.glob('/home/kmt/primatediversity/data/variants/*.snps.vcf.gz')

# chromosome names:
chromosomes = [line.split()[0] for line in open('metadata/panu3_chrom_sizes.txt').readlines()]

# index VCF files:
tabix_indexed = gwf.map(tabix_index, vcf_files)

# merge VCF files to get one merged file per chromosome/contig
merged_vcfs = []
for chrom in chromosomes:

    # for now we skip unassigned contigs
    if chrom.startswith('chrUn'):
        continue

    os.makedirs('steps/merged_vcfs', exist_ok=True)
    output_path = f'steps/merged_vcfs/{chrom}.vcf.gz'

    task = gwf.target_from_template(
        f'merge_vcf_{chrom}',
        merge_vcf_by_chrom(
            vcf_paths=vcf_files,
            tabix_paths=collect(tabix_indexed.outputs, ['path'])['paths'],
            output_path=output_path,
            chrom=chrom
        )
    )
    merged_vcfs.append(task.outputs)

# index the merged VCF files:
merged_tabix_indexed = gwf.map(tabix_index, merged_vcfs, name='tabix_merged')


