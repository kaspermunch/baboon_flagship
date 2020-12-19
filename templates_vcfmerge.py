from gwf import AnonymousTarget

def tabix_index(path):
    """Makes a tabix index on a VCF files. Existing files are overwritten.

    Args:
        path (str): Path to VCF file.

    Returns:
        gwf.AnonymousTarget: GWF target.
    """
    inputs = {'path': path}
    outputs = {'path': path + '.tbi'}
    options = {'memory': '4g',
               'walltime': '0-01:00:00'}
    spec = f'tabix -f -p vcf {path}'
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def merge_vcf_by_chrom(vcf_paths, tabix_paths, output_path, chrom):
    """Merges variants for a chromosome from multiple VCF files into a single file.

    Args:
        vcf_paths (list): List of strings representing paths to VCF files to merge.
        tabix_paths (list): List of strings representing paths to corresponding tabix files
        output_path (str): Path to merged output file.
        chrom (str): Chromosome to extract variants for.

    Returns:
        gwf.AnonymousTarget: GWF target.
    """
    inputs = {'paths': vcf_paths + tabix_paths}
    outputs = {'path': output_path}
    options = {'memory': '15g',
               'walltime': '0-03:00:00'}
    spec = f"bcftools merge --regions {chrom} --missing-to-ref -Oz -o {output_path} {' '.join(vcf_paths)}"
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)