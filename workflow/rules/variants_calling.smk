rule haplotype_caller_gene:
    input:
        bam="marked_dedup/{sample}_ref.bam",
        index = "marked_dedup/{sample}_ref.bam.bai",
        ref=config['ref']['genome'],
        l_genes=config['ref']['list'],
    output:
        gvcf="calls/{sample}_ref.g.vcf"
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/haplotypecaller/{sample}.log"
    params:
        extra="",  # optional
        java_opts=config['java_opts_parallel'],
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    script:
        "scripts/haplotypecaller_gene.py"


rule combine_gvcfs:
    input:
        gvcfs=expand("calls/{sample}_ref.g.vcf", sample = samples),
        ref=config['ref']['genome']
    output:
        gvcf="calls/all_ref_g.vcf"
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/combinegvcfs/combinegvcfs.log"
    params:
        extra = "",
        java_opts=config['java_opts_combine'],
    threads: config['threads_combine']
    resources:
        mem_mb=config['mem_mb_combine']
    script:
        "scripts/combinegvcfs.py"


rule genotype_gvcfs:
    input:
        gvcf="calls/all_ref_g.vcf",  # combined gvcf over multiple samples
        ref=config['ref']['genome']
    output:
        vcf="calls/all_ref.vcf",
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/genotypegvcfs/genotypegvcfs.log"
    params:
        extra="",  # optional
        java_opts=config['java_opts'],
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    script:
        "scripts/genotypegvcfs.py"


rule select_filter_variants:
    input:
        vcf="calls/all_ref.vcf",
        ref=config['ref']['genome']
    output:
        vcf="calls/selected.vcf"
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/select/snvs_Nile.log"
    params:
        extra="--restrict-alleles-to BIALLELIC",  # optional filter arguments, see GATK docs
        java_opts=config['java_opts'],
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    script:
        "scripts/selectvariants.py"

