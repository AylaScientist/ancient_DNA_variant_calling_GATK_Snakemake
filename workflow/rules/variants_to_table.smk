rule gatk_variantstotable:
    input:
        vcf="annotated/annotated_all_snps_"+config['params']['annotation']['output'],
        ref=config['ref']['genome']
    output:
        vcf="variants/GT_SNP.table"
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/var2table/snvs.log"
    params:
        o1="annotated/annotated_all_snps_"+config['params']['annotation']['output'],
        extra="-SMA TRUE -F CHROM -F POS -F Gene.refGene -F Func.refGene -F ExonicFunc.refGene -F AF -GF GT",  # optional filter arguments, see GATK docs
        java_opts=config['java_opts'],
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    script:
        "scripts/variantstotable.py"

