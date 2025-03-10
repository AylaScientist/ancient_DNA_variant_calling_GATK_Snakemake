# Call SNPs using ANGSD:

rule run_angsd:
    input:
        bam=expand("marked_dedup/{sample}_ref.bam", sample = samples),
        bam_index=expand("marked_dedup/{sample}_ref.bam.bai", sample = samples),
        ref=config['ref']['genome'],
        sites=config['ref']['bed'],
    output:
        saf_idx="angsd/output.saf.idx",      # SAF index file
        sfs="angsd/output.sfs",              # Site frequency spectrum
        mafs="angsd/output.mafs",            # Major/minor allele frequencies
        glf="angsd/output.glf.gz"            # Genotype likelihoods
    params:
        min_mapq=20,                   # Minimum mapping quality
        min_baseq=15                   # Minimum base quality
    log:
        "logs/angsd_snp_calling.log"
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    script:
        "scripts/angsd.py"

    