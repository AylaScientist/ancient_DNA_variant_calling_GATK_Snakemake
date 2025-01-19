rule map_damage:
    input:
        ref="genome.fasta",
        bam="bam/{sample}.bam",
    output:
        log="results/{sample}/Runtime_log.txt",  # output folder is infered from this file, so it needs to be the same folder for all output files
        GtoA3p="results/{sample}/3pGtoA_freq.txt",
        CtoT5p="results/{sample}/5pCtoT_freq.txt",
        dnacomp="results/{sample}/dnacomp.txt",
        frag_misincorp="results/{sample}/Fragmisincorporation_plot.pdf",
        len="results/{sample}/Length_plot.pdf",
        lg_dist="results/{sample}/lgdistribution.txt",
        misincorp="results/{sample}/misincorporation.txt",
    #   rescaled_bam="results/{sample}.rescaled.bam", # uncomment if you want the rescaled BAM file
    params:
        extra="--no-stats",  # optional parameters for mapdamage2 (except -i, -r, -d, --rescale)
    log:
        "logs/{sample}/mapdamage2.log",

    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    script:
        "scripts/mapdamage2.py"


rule add_replace_rg:
    input:
        bam="bam/{sample}.bam",
    output:
        "fixed-rg/{sample}_ref.bam"
    conda:
        "envs/picard.yaml"
    log:
        "logs/picard/replace_rg/{sample}.log"
    params:
        extra = "SORT_ORDER=coordinate RGID=NextSeq RGLB=idp RGPL=illumina RGPU={sample} RGSM={sample} CREATE_INDEX=True",
        java_opts=config['java_opts_parallel'],
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    script:
        "scripts/readgroups.py"


rule mark_duplicates:
    input:
        "fixed-rg/{sample}_ref.bam"
    output:
        bam="marked_dedup/{sample}_ref.bam",
        metrics="marked_dedup/{sample}_ref.metrics.txt"
    conda:
        "envs/picard.yaml"
    log:
        "logs/picard/dedup/{sample}.log"
    params:
        extra = "REMOVE_DUPLICATES=true", #Duplicates can also be removed with UMI-tools
        java_opts=config['java_opts_parallel'],
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    script:
        "scripts/mark_duplicates.py"
