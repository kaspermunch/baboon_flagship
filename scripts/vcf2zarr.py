
import sys
import os
import re
import allel
import zarr
import numpy as np

script_name, zarr_path, *vcf_paths = sys.argv

# build zarr disk data structure for call set:

for vcf_path in vcf_paths:
    # get chromosome name from file name
    chrom = re.search(r'(chr[^.]+)', vcf_path).group(1)
    # add group for chromosome
    allel.vcf_to_zarr(vcf_path, zarr_path, group=chrom, fields='*', log=sys.stdout)

# open callability masks for reading
callability_masks = zarr.open_group('./steps/called_interv.zarr', mode='r')

# open callset for reading + writing
callset = zarr.open_group(zarr_path, mode='r+')

# When we merged the VCF files, we treated missing positions in some individuals as 
# REF/REF although they could also be uncalled positions. So we need to check, for 
# each individual, if any positions in the merged VCF file are not in the callability
# mask. Below we find such positions and change the site from ref/ref to missing data
# (./. represented as (-1, -1)):

# iterate over chromosomes
for chrom, callset_chrom in callset.groups():

    # iterate over samples
    for sample_idx, sample in enumerate(callset_chrom['samples'][:]):

        # positions of snps
        pos = allel.SortedIndex(callset_chrom['variants/POS'])

        # starts and ends of called intervals
        starts, ends = callability_masks[f'{sample}/{chrom}'][:].T

        # find snp positions that overlap callabilty mask
        is_called = np.digitize(pos-1, starts) == np.digitize(pos-1, ends) + 1  # -1 to pos to account for that VCF is one based and bed is zero-based

        # get index of the ones not in the callabilty mask
        idx_not_called = np.nonzero(~is_called)[0]
        
        # get genotypes for individual
        gt = callset_chrom['calldata/GT']
    
        # change the missing ones from REF/REF to ./. (-1, -1)
        gt.set_orthogonal_selection((idx_not_called, [sample_idx], slice(None)), -1)
