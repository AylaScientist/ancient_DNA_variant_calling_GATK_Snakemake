rule gatk_variantstotable_ref:
    input:
        vcf="annotated/annotated_all_snps_"+config['params']['annotation']['output'],
        ref=config['ref']['genome']
    output:
        vcf="variants/AD_GT_counts_bi_ref.table"
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/var2table/snvs_ref.log"
    params:
        o1="annotated/annotated_all_snps_"+config['params']['annotation']['output'],
        extra="-SMA TRUE -F CHROM -F POS -F Gene.refGene -F Func.refGene -F ExonicFunc.refGene -F AF -GF AD -GF GT",  # optional filter arguments, see GATK docs
        java_opts=config['java_opts'],
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    script:
        "scripts/variantstotable.py"

