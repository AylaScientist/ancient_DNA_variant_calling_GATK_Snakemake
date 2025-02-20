rule bwa_mem:
    input:
        reads="trimmed/{sample}.1.fastq",
        idx=config['ref']["bwa_idx"]
    output:
        "bam/{sample}.bam",
    log:
        "logs/bwa_mem/{sample}.log",
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sorting="samtools",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra="",  # Extra args for samtools/picard.
        java_opts=config['java_opts_parallel']
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    script:
        "scripts/bwa_mem.py"


rule samtools_flagstat:
    input:
        "bam/{sample}.bam",
    output:
        "bam/{sample}.bam.flagstat",
    log:
        "samtools/flagstat/{sample}.log",
    params:
        extra="",  # optional params string
    threads: 1
    resources:
        mem_mb=8
    script:
        "scripts/samtools_flagstat.py"
