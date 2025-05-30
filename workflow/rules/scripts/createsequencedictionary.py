__author__ = "Johannes Köster"
__copyright__ = "Copyright 2018, Johannes Köster"
__email__ = "johannes.koester@protonmail.com"
__license__ = "MIT"


import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

extra = snakemake.params.get("extra", "")
java_opts = snakemake.params.get("java_opts")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

with tempfile.TemporaryDirectory(dir="./tmpdir") as tmpdir:
    shell(
        "picard CreateSequenceDictionary"
        " {java_opts} {extra}"
        " --REFERENCE {snakemake.input.gen}"
        " --TMP_DIR {tmpdir}"
        " --OUTPUT {snakemake.output.dicti}"
        " {log}"
    )