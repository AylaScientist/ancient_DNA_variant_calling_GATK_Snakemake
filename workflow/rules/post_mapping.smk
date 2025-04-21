rule map_damage:
    input:
        ref=config['ref']['genome'],
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

"""
rule add_replace_rg:
    input:
        "bam/{sample}.bam",
    output:
        "fixed-rg/{sample}.bam"
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
"""


rule validate_bam:
    input:
        "bam/{sample}.bam"
    output:
        "logs/validate_out/{sample}.log"
    conda:
        "envs/picard.yaml"
    log:
        "logs/validate/{sample}.log"
    threads: 1  # Validation doesn't need multi-threading
    resources:
        mem_mb=8  # Adjust based on available memory
    shell:
        """
        set -euo pipefail

        picard ValidateSamFile \
            I={input} \
            O={output} \
            MODE=SUMMARY \
            VALIDATION_STRINGENCY=STRICT \
            2> {log}
        """


rule mark_duplicates:
    input:
        "bam/{sample}.bam"
    output:
        bam="marked_dedup/{sample}.bam",
        metrics="marked_dedup/{sample}.metrics.txt"
    conda:
        "envs/picard.yaml"
    log:
        "logs/picard/dedup/{sample}.log"
    params:
        extra = "REMOVE_DUPLICATES=false", #Duplicates can also be removed with UMI-tools
        java_opts=config['java_opts_parallel'],
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    script:
        "scripts/mark_duplicates.py"

rule bam2sam:
    input:
        bam="marked_dedup/{sample}.bam",
    output:
        sam="pmd_temp/temp{sample}.sam",
    log:
        "logs/samtools/bam2sam/{sample}.log"
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    shell:
       """
        # Convert BAM to SAM and remove reads with malformed quality scores
        samtools view -h {input.bam} | awk 'length($10) == length($11)' > {output.sam}
        """

rule pmdtools:
    input:
        sam="pmd_temp/temp{sample}.sam",
    output:
        pmd="pmd_temp/temp{sample}_filtered.sam",
    log:
        "logs/pmdtools/{sample}.log"
    params:
        options="--threshold=3 --header",  # Modify options as needed 
        #temp_sam="pmd_temp/temp{sample}.sam",
        #temp_filtered="pmd_temp/temp{sample}_filtered.sam",
        temp_fixed="pmd_temp/temp{sample}_fixed.sam",
        stats="filtered_pmd/{sample}.stats"
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    shell:
       """
        # Run PMDTools, logging errors
        python /vivianelabfs/ayla/envs/Snakemake-5.30.1/bin/pmdtools {params.options} < {input.sam} > {output.pmd}
        """


rule pmdtools_deamination:
    input:
        sam="pmd_temp/temp{sample}.sam",
    output:
        pmd="pmd_temp/temp{sample}_stats.sam",
    log:
        "logs/pmdtools/{sample}.log"
    params:
        options="--deamination",  # Modify options as needed 
        #temp_sam="pmd_temp/temp{sample}.sam",
        #temp_filtered="pmd_temp/temp{sample}_filtered.sam",
        temp_fixed="pmd_temp/temp{sample}_fixed.sam",
        stats="filtered_pmd/{sample}.stats"
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    shell:
       """
        # Run PMDTools, logging errors
        python /vivianelabfs/ayla/envs/Snakemake-5.30.1/bin/pmdtools {params.options} < {input.sam} > {output.pmd}
        """

"""
rule fix_sam:
    input:
        pmd="pmd_temp/temp{sample}_filtered.sam",
    output:
        fix="pmd_temp/temp{sample}_fixed.sam",
    log:
        "logs/awk/fixed_pmd_{sample}.log"
    params:
        stats="filtered_pmd/{sample}.stats"
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    shell:
        # Convert the generated file into a sam file and a stats table
        #awk '{{if ($0 ~ /^@/) print > "{input.pmd}"; else print > "{params.stats}"}}' {input.pmd} > {output.fix}
        
"""

rule sam2bam:
    input:
        fix="pmd_temp/temp{sample}_filtered.sam",
    output:
        pmd="filtered_pmd/{sample}.bam",
    log:
        "logs/samtools/sam2bam/{sample}.log"
    params:
        extra="",
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    shell:
        """ 
        # Convert filtered SAM back to BAM
        samtools view -bo {output.pmd} {input.fix}
        """



rule samtools_index:
    input:
        "bam/{sample}.bam",
    output:
        "bam/{sample}.bam.bai",
    log:
        "logs/samtools_index/{sample}.log",
    params:
        extra="",  # optional params string
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    script:
        "scripts/samtools_index.py"
