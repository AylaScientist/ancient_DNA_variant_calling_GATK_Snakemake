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


rule genome_dictionary:
    input:
        gen = config['ref']['genome'],
    output:
        dicti = config['ref']["dict"],
    log:
        "logs/picard/create_seq_dictionary.log",
    params:
        java_opts=config['java_opts'],
        extra="", # optional: extra arguments for picard.
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    script:
        "scripts/createsequencedictionary.py"


rule samtools_fai:
    input:
        config['ref']['genome'],
    output:
        config['ref']['fai'],
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    shell:
        (f"samtools faidx {input}")


rule fastp_se:
    input:
        sample=["fastq_merged/{sample}.fastq.gz"]
    output:
        trimmed="trimmed/se/{sample}.fastq",
        failed="trimmed/se/{sample}.failed.fastq",
        html="trimmed/report/se/{sample}.html",
        json="trimmed/report/se/{sample}.json"
    log:
        "logs/fastp/se/{sample}.log"
    params:
        adapters="", #"--adapter_sequence ACGGCTAGCTA",
        extra=""
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    script:
        "scripts/fastp.py"


rule trimmomatic_filter:
    input:
        r1 = "trimmed/se/{sample}.fastq",
    output:
        r1 = "trimmed/{sample}.1.fastq",
        # reads where trimming entirely removed the mate
        #r1_unpaired = temp("trimmed/{sample}.1.unpaired.fastq"),
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
        tempdir = "./trimmomatic_tempdir",
        r1_unpaired = "trimmed/{sample}.1.unpaired.fastq",
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    script:
        "scripts/trimmomatic_SE.py"
