__author__ = "Aurora Campo"
__copyright__ = "Copyright 2024, Aurora Campo"
__license__ = "MIT"

import os.path
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

in_bam = snakemake.input.get("bam", "")
if in_bam:
    in_bam = "--input " + in_bam


shell:(
    "angsd -doSaf 1"
    " -doMajorMinor 1"
    " -doCounts 1"
    " -GL 2 "
    " -bam {snakemake.input.bam}"
    " -anc {snakemake.input.ref}"
    " -out output "
    " -minMapQ {snakemake.params.min_mapq}"
    " -minQ {snakemake.params.min_baseq}"
    " -sites {snakemake.input.sites}"
    " -nThreads {snakemake.threads}"
    " &> {log}"
    "realSFS {snakemake.output.saf<.idx} > {snakemake.output.sfs}"
)