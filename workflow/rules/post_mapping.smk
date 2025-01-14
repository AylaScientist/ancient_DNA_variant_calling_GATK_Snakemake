rule map_damage:
    input:
        reads="processed_reads/trimmed.fq"
    output:
        report="damage_reports/report_{wildcards.sample}.txt"
    shell:
        """
        mapDamage -i {input.reads} -r {REF}
        """

rule add_replace_rg:
    input:
        bam="mapped_reads/aligned_{wildcards.sample}.bam"
    output:
        rg_bam="mapped_reads/rg_aligned_{wildcards.sample}.bam"
    shell:
        """
        picard AddOrReplaceReadGroups I={input.bam} O={output.rg_bam} \
        RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM={wildcards.sample}
        """

rule mark_duplicates:
    input:
        bam="mapped_reads/rg_aligned_{wildcards.sample}.bam"
    output:
        dedup_bam="processed_reads/dedup_{wildcards.sample}.bam",
        metrics="processed_reads/{wildcards.sample}_metrics.txt"
    shell:
        """
        picard MarkDuplicates I={input.bam} O={output.dedup_bam} \
        M={output.metrics} REMOVE_DUPLICATES=true
        """
