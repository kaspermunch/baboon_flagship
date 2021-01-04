import sys, os
import numpy as np
from collections import defaultdict
import zarr

script_name, callable_zarr_path, callset_zarr_path, chromosome_lengths_path, callmasks_zarr_path = sys.argv

# called intervals zarr
callability_masks = zarr.open_group(callable_zarr_path, mode='r') 

# callset zarr
callset = zarr.open_group(callset_zarr_path, mode='r') 

# chromosome lengths
chromosome_lengths = dict()
for line in open(chromosome_lengths_path):
    chrom, length = line.split()
    chromosome_lengths[chrom] = int(length)

# output path for callmask zarr
root = zarr.open_group(callmasks_zarr_path, mode='w')

# loop chromosomes
for chrom, length in chromosome_lengths.items():

    # get samples for chromosome
    samples = callset[f'{chrom}/samples'][:]

    # add big matrix for all individuals
    callability = root.array(chrom, np.zeros((len(samples), length), dtype='b'))

    # loop samples
    for sample_idx, sample in enumerate(samples):

        # get called positions for individual
        called_positions = []
        for start, end in callability_masks[f'{sample}/{chrom}'][:]:
            called_positions.extend(range(start, end))

        # update matrix for all individuals with called positions
        callability.oindex[[sample_idx], called_positions] = 1



# # loop chromosomes
# for chrom, length in chromosome_lengths.items():

#     # get samples for chromosome
#     samples = callset[f'{chrom}/samples'][:]

#     # add group for chromosome
#     group = root.create_group(chrom)

#     # add big matrix for all individuals
#     callability_all = group.array('ALL', np.zeros((len(samples), length), dtype='b'))

#     # loop samples
#     for sample_idx, sample in enumerate(samples):

#         # get called positions for individual
#         called_positions = []
#         for start, end in callability_masks[f'{sample}/{chrom}'][:]:
#             called_positions.extend(range(start, end))

#         # add array for individual
#         callability = group.array(sample, np.zeros(length), dtype='b')

#         # update array for individual with called positions
#         callability.set_orthogonal_selection(called_positions, 1)

#         # update matrix for all individuals with called positions
#         callability_all.set_orthogonal_selection(([sample_idx], called_positions), 1)            


