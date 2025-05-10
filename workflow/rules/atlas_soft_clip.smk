rule filter_bam:
    input:
        "atlas/{sample}.trimmed_softClippedBasesRemoved.bam"
    output:
        "atlas/{sample}_filtered.bam"
    params:
        mapq=config["params"]["atlas"]["min_mapq"],
        length=config["params"]["atlas"]["min_read_length"],
        prefix='atlas/{sample}'
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    log: "logs/atlas/filter_bam/{sample}.log"
    shell:
        """
        atlas task=filter bam=./{input} minMapQ={params.mapq} minReadLength={params.length} out={params.prefix} > {log} 2>&1
        """

rule extract_region:
    input:
        bam="atlas/{sample}_filtered.bam"
    output:
        "atlas/{sample}_region.bam"
    params:
        region=config["params"]["atlas"]["region"] #area to call SNPs from
    shell:
        """
        samtools view -b {input.bam} {params.region} > {output}
        samtools index {output}
        """


rule allelic_depth:
    input:
        "atlas/{sample}_region.bam"
    output:
        "atlas/{sample}_allelicDepth.txt"
    params:
        prefix="atlas/{sample}",
        depth=[config][params][atlas][depth]
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    log: "logs/atlas/allelic_depth/{sample}.log"
    shell:
        """
        atlas task=allelicDepth bam={input} readUpToDepth={params.depth} out={params.prefix} > {log} 2>&1
        """

rule fix_bug_atlas:
    input:
        "atlas/{sample}_allelicDepth.txt"
    output:
        "atlas/{sample}_allelicDepth.txt.gz"
    shell:
        """
        mv {input} {output}
        """

rule call_snps:
    input:
        glf="atlas/{sample}_allelicDepth.txt.gz",
        bam="atlas/{sample}_region.bam",
        ref=config['ref']['genome']
    output:
        "calls_atlas/{sample}_MaximumLikelihood.vcf.gz"
    params:
        prefix="calls_atlas/{sample}"
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    log: "logs/atlas/call_snps/{sample}.log"
    shell:
        """
        atlas task=call bam={input.bam} glf={input.glf} fasta={input.ref} out={params.prefix} > {log} 2>&1
        """