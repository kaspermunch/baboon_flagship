
import sys
import os
import re
import allel
import zarr
from collections import defaultdict
import gzip

script_name, zarr_path, *bed_paths = sys.argv

root = zarr.open_group(zarr_path, mode='w')

for bed_file in bed_paths:
    sample = os.path.basename(bed_file).split('.')[0]

    coordinates = defaultdict(list)
    with gzip.open(bed_file, 'rt') as f:
        for line in f:
            chrom, start, end = line.split()
            coordinates[chrom].append((int(start), int(end)))

    indiv = root.create_group(sample)
    for chrom in sorted(coordinates):
        indiv.array(chrom, coordinates[chrom], dtype='int32')

