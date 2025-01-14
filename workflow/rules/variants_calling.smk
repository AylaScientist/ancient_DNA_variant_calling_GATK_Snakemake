rule haplotype_caller:
    input:
        bam="processed_reads/dedup_{wildcards.sample}.bam"
    output:
        gvcf="variants/{wildcards.sample}.g.vcf"
    shell:
        """
        gatk --java-options "-Xmx4g" HaplotypeCaller \
        -R {REF} -I {input.bam} -O {output.gvcf} -ERC GVCF
        """

rule genotype_gvcfs:
    input:
        gvcfs=expand("variants/{{sample}}.g.vcf", sample=SAMPLES)
    output:
        vcf="variants/combined.vcf"
    shell:
        """
        gatk --java-options "-Xmx4g" GenotypeGVCFs \
        -R {REF} -V {input.gvcfs} -O {output.vcf}
        """

rule select_filter_variants:
    input:
        vcf="variants/combined.vcf"
    output:
        filtered_vcf="variants/filtered_combined.vcf"
    shell:
        """
        gatk SelectVariants \
        -R {REF} -V {input.vcf} -select-type SNP \
        -O {output.filtered_vcf}
        """
