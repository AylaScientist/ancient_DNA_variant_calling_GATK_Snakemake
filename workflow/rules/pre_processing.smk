rule bwa_index:
    input:
        genome = config['ref']['genome'],
    output:
        idx = config['ref']["bwa_idx"],
    log:
        "logs/bwa_index/bwa_genome_idex.log",
    params:
        extra=lambda w: f"-a bwtsw",
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    script:
        "scripts/bwa_index.py"


rule adapter_removal_se:
    input:
        sample = ["fastq_merged/{sample}.fastq.gz"],
    output:
        fq="trimmed/se/{sample}.fastq.gz",                               # trimmed reads
        discarded="trimmed/se/{sample}.discarded.fastq.gz",              # reads that did not pass filters
        settings="adapter_stats/se/{sample}.settings"                    # parameters as well as overall statistics
    log:
        "logs/adapterremoval/se/{sample}.log"
    params:
        adapters=config['params']['adapters'],
        extra="",
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    script:
        "scripts/adapter_removal_SE.py"


rule trimmomatic_filter:
    input:
        r1 = "trimmed/se/{sample}.fastq.gz",
    output:
        r1 = ("trimmed/{sample}.1.fastq"),
        # reads where trimming entirely removed the mate
        r1_unpaired = temp("trimmed/{sample}.1.unpaired.fastq"),
    conda:
        "envs/trimmomatic.yaml"
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer = config['params']['trimmomatic']['se']['trimmer'], # optional parameters
        #"ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:5 LEADING:10 TRAILING:10 MINLEN:38",
        extra = "",
        java_opts = "", # config['java_opts_parallel'],
        path = config['params']['trimmomatic']['path'],
        tempdir = "./trimmomatic_tempdir"
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    script:
        "scripts/trimmomatic_SE.py"