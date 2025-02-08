
rule bwa_aln:
    input:
        fastq="trimmed/{sample}.1.fastq",
        idx=config['ref']["bwa_idx"],
    output:
        "sai/{sample}.sai",
    params:
        extra="",
    log:
        "logs/bwa_aln/{sample}.log",
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    script:
        "scripts/bwa_aln.py"


rule bwa_samse:
    input:
        fastq="trimmed/se/{sample}.fastq.gz",
        sai="sai/{sample}.sai",
        idx=config['ref']["bwa_idx"],
    output:
        "bam/{sample}.bam",
    params:
        extra=r"-r '@RG\tID:{sample}\tSM:{sample}'",  # optional: Extra parameters for bwa.
        sort="samtools",  # optional: Enable sorting. Possible values: 'none', 'samtools' or 'picard'`
        sort_order="coordinate",  # optional: Sort by 'queryname' or 'coordinate'
        sort_extra="",  # optional: extra arguments for samtools/picard
    log:
        "logs/bwa_samse/{sample}.log",
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    script:
        "scripts/bwa_samse.py"