"""
bio/trimmomatic/pe
Snakemake wrapper to trim reads with trimmomatic in PE mode with help of pigz.
pigz is the parallel implementation of gz. Trimmomatic spends most of the time
compressing and decompressing instead of trimming sequences. By using process
substitution (<(command), >(command)), we can accelerate trimmomatic a lot.
Consider providing this wrapper with at least 1 extra thread per each gzipped
input or output file.
"""

__author__ = "Aurora Campo modified from Johannes Köster, Jorge Langa"
__copyright__ = "Copyright 2025, Johannes Köster"
__email__ = "ayla.bcn@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

# Distribute available threads between trimmomatic itself and any potential pigz instances
def distribute_threads(input_files, output_files, available_threads):
    gzipped_input_files = sum(1 for file in input_files if file.endswith(".gz"))
    gzipped_output_files = sum(1 for file in output_files if file.endswith(".gz"))
    potential_threads_per_process = available_threads // (
        1 + gzipped_input_files + gzipped_output_files
    )
    if potential_threads_per_process > 0:
        # decompressing pigz creates at most 4 threads
        pigz_input_threads = (
            min(4, potential_threads_per_process) if gzipped_input_files != 0 else 0
        )
        pigz_output_threads = (
            (available_threads - pigz_input_threads * gzipped_input_files)
            // (1 + gzipped_output_files)
            if gzipped_output_files != 0
            else 0
        )
        trimmomatic_threads = (
            available_threads
            - pigz_input_threads * gzipped_input_files
            - pigz_output_threads * gzipped_output_files
        )
    else:
        # not enough threads for pigz
        pigz_input_threads = 0
        pigz_output_threads = 0
        trimmomatic_threads = available_threads
    return trimmomatic_threads, pigz_input_threads, pigz_output_threads


def compose_input_gz(filename, threads):
    if filename.endswith(".gz") and threads > 0:
        return "<(pigz -p {threads} --decompress --stdout {filename})".format(
            threads=threads, filename=filename
        )
    return filename


def compose_output_gz(filename, threads, compression_level):
    if filename.endswith(".gz") and threads > 0:
        return ">(pigz -p {threads} {compression_level} > {filename})".format(
            threads=threads, compression_level=compression_level, filename=filename
        )
    return filename

def get_java_opts(snakemake):
    java_opts = ""
    mem_mb = snakemake.resources.get("mem_mb", 4000) - snakemake.params.get("java_mem_overhead_mb", 256)
    java_opts += f" -Xmx{mem_mb}M"
    tmpdir = snakemake.resources.get("tmpdir", "/tmp")  # Use a default /tmp directory if not provided
    java_opts += f" -Djava.io.tmpdir={tmpdir}"
    return java_opts

extra = snakemake.params.get("extra", "")
java_opts = get_java_opts(snakemake)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
compression_level = snakemake.params.get("compression_level", "-5")
trimmer = " ".join(snakemake.params.trimmer)

# Distribute threads
input_files = [snakemake.input.r1]
output_files = [
    snakemake.output.r1,
]

trimmomatic_threads, input_threads, output_threads = distribute_threads(
    input_files, output_files, snakemake.threads
)

input_r1 = [
    compose_input_gz(filename, input_threads) for filename in input_files
]

output_r1 = [
    compose_output_gz(filename, output_threads, compression_level)
    for filename in output_files
]

shell(
    "java {java_opts} -jar /powerapps/share/centos7/trimmomatic/0.39/trimmomatic-0.39.jar SE "
    "-threads {trimmomatic_threads} -phred33 {extra} "
    "{snakemake.input.r1} "
    "{snakemake.output.r1} "# {output_r1_unp} "
    "{trimmer} "
    "{log}"
)
