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
        cp {log} {output}
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


rule pmdtools:
    input:
        bam="marked_dedup/{sample}.bam",
    output:
        filtered_bam="filtered_pmd/{sample}.bam",
    log:
        "logs/pmdtools/{sample}.log"
    params:
        options="--deamination --threshold=3 --header",  # Modify options as needed 
        temp_sam="pmd_temp/temp{sample}.sam",
        temp_filtered="pmd_temp/temp{sample}_filtered.sam"
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    shell:
       """
        # Convert BAM to SAM and remove reads with malformed quality scores
        samtools view -h {input.bam} | awk 'length($10) == length($11)' > > {params.temp_sam}

        # Run PMDTools, logging errors
        python /vivianelabfs/ayla/envs/Snakemake-5.30.1/bin/pmdtools {params.options} < {params.temp_sam} > {params.temp_filtered} 2> {log}

        # Ensure PMDTools produced output before proceeding
        if [ ! -s {params.temp_filtered} ]; then
            echo "Error: PMDTools did not generate a valid SAM file" >> {log}
            exit 1
        fi

        # Convert filtered SAM back to BAM
        samtools view -bo {output.filtered_bam} {params.temp_filtered}

        # Clean up
        rm temp.sam {params.temp_filtered}
        """


rule samtools_index:
    input:
        "filtered_pmd/{sample}.bam",
    output:
        "filtered_pmd/{sample}.bam.bai",
    log:
        "logs/samtools_index/{sample}.log",
    params:
        extra="",  # optional params string
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    script:
        "scripts/samtools_index.py"
