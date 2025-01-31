__author__ = "Aurora Campo"
__copyright__ = "Copyright 2025, Aurora Campo"
__license__ = "MIT"

import os.path
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Ensure output directory exists
out_prefix = os.path.dirname(snakemake.output.saf_idx)
os.makedirs(out_prefix, exist_ok=True)

# Create a properly formatted BAM list file
bam_list_file = os.path.join(out_prefix, "bamlist.txt")

with open(bam_list_file, "w") as f:
    for bam in snakemake.input.bam:
        f.write(f"{bam}\n")  # Ensure one BAM per line

# Print BAM list for debugging
print(f"\n[INFO] BAM list written to {bam_list_file}")
with open(bam_list_file, "r") as f:
    print(f.read())


# Ensure the BED file is indexed
bed_file = snakemake.input.sites
bed_index_file = f"{bed_file}.bin"

if not os.path.exists(bed_index_file):
    print(f"\n[INFO] Indexing BED file: {bed_file}")
    shell(f"angsd sites index {bed_file}")
else:
    print(f"\n[INFO] BED file index already exists: {bed_index_file}")


# Check and fix FASTA index timestamp issue
fasta_file = snakemake.input.ref
fasta_index = f"{fasta_file}.fai"

if not os.path.exists(fasta_index) or os.path.getmtime(fasta_index) < os.path.getmtime(fasta_file):
    print(f"\n[INFO] Re-indexing FASTA file: {fasta_file}")
    shell(f"samtools faidx {fasta_file}")
else:
    print(f"\n[INFO] FASTA index is up-to-date: {fasta_index}")


"""
ANGSD supports different genotype likelihood models:
-GL     Value	Model-Used-Notes
-GL     1       SAMtools	Standard for BAM files (works well for ancient DNA).
-GL     2	    GATK	Uses GATK's genotype likelihood model.
-GL     3	    SOAPsnp	Rarely used, mostly for short reads.
-GL     4	    SYK model	Sanger sequencing-specific.

💡 Recommendation: Use -GL 2 (GATK model) for ancient DNA unless you have a reason to use -GL 1.
"""


# Construct the ANGSD command
angsd_cmd = (
    f"angsd "
    f" -doGlf 2 -GL 2"  # Output genotype likelihoods in Beagle format., Chosen GL2 as GATK likelihood model
    f" -doMajorMinor 1"
    f" -doMaf 2"  # Output individual-level minor allele frequencies.
    f" -doPost 1"  # Calculate posterior genotype probabilities from likelihoods.
    f" -bam {bam_list_file}"
    f" -anc {snakemake.input.ref}"
    f" -out {snakemake.output.saf_idx.rsplit('.', 1)[0]}"  # Match ANGSD output prefix
    f" -minMapQ {snakemake.params.min_mapq}"
    f" -minQ {snakemake.params.min_baseq}"
    f" -sites {snakemake.input.sites}"
    f" -nThreads {snakemake.threads} "
    f" {log}"
)

# Print the command for logging/debugging
print("\n[INFO] Running ANGSD with the following command:\n")
print(angsd_cmd + "\n")

# Run the ANGSD command
shell(angsd_cmd)

# Construct and print the realSFS command
realSFS_cmd = f"realSFS {snakemake.output.saf_idx} > {snakemake.output.sfs}"

print("\n[INFO] Running realSFS with the following command:\n")
print(realSFS_cmd + "\n")

# Run realSFS to generate the site frequency spectrum
shell(realSFS_cmd)