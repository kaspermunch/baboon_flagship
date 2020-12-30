import os, sys
from gwf import Workflow
from gwf.workflow import collect
from templates_vcfmerge import *

gwf = Workflow(defaults={'account': 'primatediversity'})


###############################################################################
# Input data
###############################################################################

# Samples with strange contig names...:
# PD_0265
# PD_0266
# PD_0267
# PD_0268
# PD_0269
# PD_0270
# PD_0271
# PD_0272

# paths to VCF files:
baboon_samples = [
 'PD_0199', 'PD_0200', 'PD_0201', 'PD_0202', 'PD_0203', 'PD_0204', 'PD_0205', 'PD_0206', 'PD_0207', 
 'PD_0208', 'PD_0209', 'PD_0210', 'PD_0211', 'PD_0212', 'PD_0213', 'PD_0214', 'PD_0215', 'PD_0216', 
 'PD_0217', 'PD_0218', 'PD_0219', 'PD_0220', 'PD_0221', 'PD_0222', 'PD_0223', 'PD_0224', 'PD_0225', 
 'PD_0226', 'PD_0227', 'PD_0228', 'PD_0229', 'PD_0230', 'PD_0231', 'PD_0232', 'PD_0233', 'PD_0234', 
 'PD_0235', 'PD_0236', 'PD_0237', 'PD_0238', 'PD_0239', 'PD_0240', 'PD_0241', 'PD_0242', 'PD_0243', 
 'PD_0244', #'PD_0265', 'PD_0266', 'PD_0267', 'PD_0268', 'PD_0269', 'PD_0270', 'PD_0271', 'PD_0272', 
 'PD_0390', 'PD_0391', 'PD_0392', 'PD_0393', 'PD_0394', 'PD_0395', 'PD_0396', 'PD_0397', 'PD_0398', 
 'PD_0399', 'PD_0400', 'PD_0492', 'PD_0401', 'PD_0493', 'PD_0494', 'PD_0495', 'PD_0496', 'PD_0497', 
 'PD_0498', 'PD_0499', 'PD_0500', 'PD_0501', 'PD_0502', 'PD_0503', 'PD_0504', 'PD_0505', 'PD_0506', 
 'PD_0507', 'PD_0508', 'PD_0509', 'PD_0637', 'PD_0641', 'PD_0642', 'PD_0649', 'PD_0650', 'PD_0651',
 'PD_0652', 'PD_0653', 'PD_0654', 'PD_0658', 'PD_0659', 'PD_0662', 'PD_0674', 'PD_0675', 'PD_0676', 
 'PD_0677', 'PD_0678', 'PD_0679', 'PD_0680', 'PD_0681', 'PD_0682', 'PD_0683', 'PD_0684', 'PD_0685', 
 'PD_0686', 'PD_0687', 'PD_0688', 'PD_0689', 'PD_0690', 'PD_0691', 'PD_0692', 'PD_0693', 'PD_0694', 
 'PD_0695', 'PD_0696', 'PD_0697', 'PD_0698', 'PD_0699', 'PD_0700', 'PD_0701', 'PD_0702', 'PD_0703', 
 'PD_0704', 'PD_0705', 'PD_0706', 'PD_0707', 'PD_0708', 'PD_0709', 'PD_0710', 'PD_0711', 'PD_0712', 
 'PD_0713', 'PD_0714', 'PD_0715', 'PD_0716', 'PD_0717', 'PD_0718', 'PD_0719', 'PD_0720', 'PD_0721', 
 'PD_0722', 'PD_0723', 'PD_0724', 'PD_0725', 'PD_0726', 'PD_0727', 'PD_0728', 'PD_0729', 'PD_0730', 
 'PD_0731', 'PD_0732', 'PD_0733', 'PD_0734', 'PD_0735', 'PD_0736', 'PD_0737', 'PD_0738', 'PD_0739', 
 'PD_0740', 'PD_0741', 'PD_0742', 'PD_0743', 'PD_0744', 'PD_0745', 'PD_0746', 'PD_0747', 'PD_0748', 
 'PD_0749', 'PD_0750', 'PD_0751', 'PD_0752', 'PD_0753', 'PD_0754', 'PD_0755', 'PD_0756', 'PD_0757', 
 'PD_0758', 'PD_0759', 'PD_0760', 'PD_0761', 'PD_0762', 'PD_0763', 'PD_0764', 'PD_0765', 'PD_0766', 
 'PD_0767', 'PD_0768', 'PD_0769', 'PD_0770', 'PD_0771', 'PD_0772', 'PD_0773', 'PD_0774', 'PD_0775', 
 'PD_0776', 'PD_0777', 'PD_0778', 'PD_0779', 'PD_0780', 'PD_0781', 'PD_0782', 'PD_0783', 'PD_0784', 
 'PD_0785', 'PD_0786', 'PD_0787', 'PD_0788', 'PD_0789', 'PD_0790', 'PD_0791', 'PD_0792', 'PD_0793',
 #'PD_0794_BAB'
 ]

vcf_files = [f'/home/kmt/primatediversity/data/variants/{sample}.variable.filtered.HF.snps.vcf.gz' for sample in baboon_samples]

bed_files = [f'/home/kmt/primatediversity/data/callability/{sample}.variable.filtered.callable.bed.gz' for sample in baboon_samples]

# lengths of chromosomes that snps are mapped to:
chromosome_lengths_path = './metadata/macFas5.chrom.sizes.txt'


## HACK FOR NOW: only include samples where both vcf and bed files exist: ----------------------
vcf_files, bed_files = zip(*[(vcf, bed) for vcf, bed in zip(vcf_files, bed_files) if os.path.exists(vcf) and os.path.exists(bed)])
vcf_files = list(vcf_files)
bed_files = list(bed_files)
assert len(vcf_files) == len(bed_files)
## ----------------------------------------------------------------------------------------------

# chromosome names:
#chromosomes = [line.split()[0] for line in open('metadata/panu3_chrom_sizes.txt').readlines()]
chromosomes = [f'chr{chrom}' for chrom in range(1,21)] + ['chrX']

###############################################################################
# Index VCF files
###############################################################################

tabix_indexed = gwf.map(tabix_index, vcf_files)

###############################################################################
# Merge VCF files to get one merged file per chromosome/contig
###############################################################################

merged_vcfs = []
for chrom in chromosomes:

    os.makedirs('steps/merged_vcfs', exist_ok=True)

    task = gwf.target_from_template(
        f'merge_vcf_{chrom}',
        merge_vcf_by_chrom(
            vcf_paths=vcf_files,
            tabix_paths=collect(tabix_indexed.outputs, ['path'])['paths'],
            output_path=f'steps/merged_vcfs/{chrom}.vcf.gz',
            chrom=chrom
        )
    )
    merged_vcfs.append(task.outputs)

###############################################################################
# Index the merged VCF files:
###############################################################################

merged_tabix_indexed = gwf.map(tabix_index, merged_vcfs, name='tabix_merged')

###############################################################################
# transform bed files to zarr: (NB: the bed2zarr.py script skips unassigned chromosomes 'chrUn')
###############################################################################

callable_zarr_path = './steps/called_interv.zarr'

gwf.target(name='bed2zarr',
     inputs=bed_files, 
     outputs=[callable_zarr_path], 
     walltime='1-00:00:00', 
     memory='8g') << f"""
mkdir -p {os.path.dirname(callable_zarr_path)}
rm -rf {callable_zarr_path}
python ./scripts/bed2zarr.py {callable_zarr_path} {' '.join(bed_files)}
"""

###############################################################################
# transform vcf into zarr:
###############################################################################

callset_zarr_path = './steps/callset.zarr'
merged_vcf_paths = collect(merged_vcfs, ['path'])['paths']

gwf.target(name='vcf2zarr',
     inputs=merged_vcf_paths + [callable_zarr_path], 
     outputs=[callset_zarr_path], 
     walltime='1-00:00:00', 
     memory='8g') << f"""
mkdir -p {os.path.dirname(callset_zarr_path)}
rm -rf {callset_zarr_path}
python ./scripts/vcf2zarr.py {callset_zarr_path} {' '.join(merged_vcf_paths)}
"""

###############################################################################
# make call masks for each base in each individual:
###############################################################################

callmasks_zarr_path = './steps/call_masks.zarr'

gwf.target(name='callmasks',
     inputs=[callset_zarr_path, callable_zarr_path, chromosome_lengths_path], 
     outputs=[callmasks_zarr_path], 
     walltime='2-00:00:00', 
     memory='24g') << f"""
mkdir -p {os.path.dirname(callmasks_zarr_path)}
rm -rf {callmasks_zarr_path}
python ./scripts/callmasks.py {callable_zarr_path} {callset_zarr_path} {chromosome_lengths_path} {callmasks_zarr_path}
"""


