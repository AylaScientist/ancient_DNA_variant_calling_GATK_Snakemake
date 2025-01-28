__author__ = "Aurora Campo"
__copyright__ = "Copyright 2024, Aurora Campo"
__license__ = "MIT"

import os.path
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell:(
    "angsd"
    " -doGlf 2" #Output genotype likelihoods in Beagle format.
    " -doMajorMinor 1" #
    " -doMaf 2" #Output individual-level minor allele frequencies.
    " -doPost 1" #Calculate posterior genotype probabilities from likelihoods.
    " -bam {snakemake.input.bam}"
    " -anc {snakemake.input.ref}"
    " -out output "
    " -minMapQ {snakemake.params.min_mapq}"
    " -minQ {snakemake.params.min_baseq}"
    " -sites {snakemake.input.sites}"
    " -nThreads {snakemake.threads}"
    " {log}"
    #"realSFS {snakemake.output.saf<.idx} > {snakemake.output.sfs}"
)