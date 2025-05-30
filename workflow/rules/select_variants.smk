rule gzip_on_vcfs:
    input:
        "calls_atlas/{sample}_MaximumLikelihood.vcf.gz"
    output:
        "calls_atlas/{sample}_MaximumLikelihood.vcf"
    params:
        backup = "calls_backup/{sample}_MaximumLikelihood.vcf.gz"
    shell:
        """
        #cp {input} {params.backup}
        gunzip {input}
        """

rule copy_vcfs:
    input:
        "calls_atlas/{sample}_MaximumLikelihood.vcf"
    output:
        "calls_atlas_bgzip/{sample}_MaximumLikelihood.vcf"
    params:
        "calls_atlas_bgzip/"
    shell:
        """
        cp {input} {params}
        """


rule bgzip_on_vcfs:
    input:
        "calls_atlas_bgzip/{sample}_MaximumLikelihood.vcf"
    output:
        "calls_atlas_bgzip/{sample}_MaximumLikelihood.vcf.gz"
    log: "logs/bgzip/{sample}.log"
    shell:
        """
        bgzip {input} > {log} 2>&1
        tabix {output}
        """    

rule merge_atlas_vcfs:
    input:
        expand("calls_atlas_bgzip/{sample}_MaximumLikelihood.vcf.gz", sample = samples),
    output:
        "calls_atlas/all_atlas.vcf.gz"
    conda:
        "envs/bcftools.yaml"
    log:
        "logs/bcftools/merge_atlas.log"
    threads: config["threads_combine"]
    resources:
        mem_mb=config["mem_mb_combine"]
    shell:
        """
        bcftools merge {input} -O z -o {output} > {log} 2>&1
        """


rule index_merged_vcf:
    input:
        "calls_atlas/all_atlas.vcf.gz"
    output:
        "calls_atlas/all_atlas.vcf.gz.tbi"
    conda:
        "envs/bcftools.yaml"
    log:
        "logs/bcftools/index_atlas.log"
    threads: config["threads_combine"]
    resources:
        mem_mb=config["mem_mb_combine"]
    shell:
        """
        bcftools index --tbi {input} > {log} 2>&1
        """


rule select_filter_variants:
    input:
        vcf="calls_atlas/all_atlas.vcf.gz",
        ref=config['ref']['genome'],
        idx="calls_atlas/all_atlas.vcf.gz.tbi"
    output:
        vcf="calls_atlas/selected_atlas.vcf.gz"
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


rule gatk_variantstotable:
    input:
        vcf="calls_atlas/selected_atlas.vcf.gz",
        ref=config['ref']['genome']
    output:
        vcf="variants_atlas/GT_SNP.table"
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/var2table/snvs.log"
    params:
        extra="-SMA TRUE -F CHROM -F POS -F AF -GF GT -GF AD",  # optional fields and flags
        java_opts=config['java_opts'],
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    shell:
        """
        gatk --java-options "{params.java_opts}" VariantsToTable \
            -V {input.vcf} \
            -R {input.ref} \
            {params.extra} \
            -O {output.vcf} > {log} 2>&1
        """

