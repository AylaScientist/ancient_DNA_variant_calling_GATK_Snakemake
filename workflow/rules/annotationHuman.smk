# Convert th vcf file to annovar format thus extracting the SNPs
rule convert_to_annovar:
    input:
        gvcf1="calls/selected.vcf"
    output:
        o1= config['params']['annotation']['convert']
    resources:
        mem_mb=config['mem_mb']
    threads: config['threads']
    conda:
        "envs/perl.yaml"
    shell:
        "perl ./rules/scripts/convert2annovar.pl {input.gvcf1} -format vcf4 -allsample -withfreq -withfilter -context -out {output.o1}"



#Make tokens for the annotation script

rule token_annotation:
    input:
        config['params']['annotation']['convert']
    output:
        o1 = config['params']['annotation']['output_annotate']
    shell:
        "touch {output.o1} "


rule token_pathbuild:
    input:
        config['params']['annotation']['convert']
    output:
        o1 = config['params']['annotation']['pathbuild']
    shell:
        "touch {output.o1} "


# Collect the annotations in the db
rule annotate_download_db:
    input:
        i1=config['params']['annotation']['convert'],
        buildver= config['params']['annotation']['buildver'],
    output:
        dbtype="humandb/" + config['params']['annotation']['buildver'] + "_" + config['params']['annotation']['dbtype']+".txt",
    params:
        i2=config['params']['annotation']['output_annotate'],
        dbtype=config['params']['annotation']['dbtype'],
        path = config['params']['annotation']['path']
    log:
        log=config['params']['annotation']['output_annotate']+".log"
    resources:
        mem_mb=config['mem_mb']
    threads: config['threads']
    conda:
        "envs/perl.yaml"
    shell:
        #"perl ./rules/scripts/annotate_variation.pl -geneanno {input.i1} -buildver {input.buildver} ./ -outfile {params.i2}" #for annotating genes
        #"perl ./rules/scripts/annotate_variation.pl -geneanno {input.i1} -buildver hg19 -downdb -webfrom annovar refGene {params.path}/ -outfile {params.i2}" #To download the db for annotating genes
        "perl ./rules/scripts/annotate_variation.pl -buildver {input.buildver} -downdb -webfrom annovar {params.dbtype} humandb/ -outfile {params.dbtype}" #To download the SNPdb 
        


# Collect the annotations in the db
rule annotate:
    input:
        i1=config['params']['annotation']['convert'],
        buildver= config['params']['annotation']['buildver'],
        dbtype="humandb/" + config['params']['annotation']['buildver'] + "_" + config['params']['annotation']['dbtype']+".txt",
    output:
        o2=config['params']['annotation']['convert'] + ".hg38_avsnp151"
    params:
        i2=config['params']['annotation']['output_annotate'],
        dbtype=config['params']['annotation']['dbtype'],
        path = config['params']['annotation']['path']
    log:
        log=config['params']['annotation']['output_annotate']+".log"
    resources:
        mem_mb=config['mem_mb']
    threads: config['threads']
    conda:
        "envs/perl.yaml"
    shell:
        #"perl ./rules/scripts/annotate_variation.pl -geneanno {input.i1} -buildver {input.buildver} ./ -outfile {params.i2}" #for annotating genes
        #"perl ./rules/scripts/annotate_variation.pl -geneanno {input.i1} -buildver hg19 -downdb -webfrom annovar refGene {params.path}/ -outfile {params.i2}" #To download the db for annotating genes
        #"perl ./rules/scripts/annotate_variation.pl -filter -buildver {input.buildver} -dbtype {params.dbtype} -out {output.o2} {input.i1} {params.path}/" #for annotating SNPs
        "perl ./rules/scripts/annotate_variation.pl -regionanno -buildver {input.buildver} -dbtype {params.dbtype} {output.o2} {input.i1} {params.path}/" #for annotating SNPs


# Make tokens for the tables with the annotation:
rule token_table:
    input:
        config['params']['annotation']['convert'] + ".hg38_avsnp151"
    output:
        "annotated/annotated_all_snps_"
    shell:
        "touch {output}"


# Annotate the file. This file can be used for the creation of the pseudogenomes
rule an_table:
    input:
        i1=config['params']['annotation']['convert'],
        i2=config['params']['annotation']['convert'] + ".hg38_avsnp151",
        i3="annotated/annotated_all_snps_",
        gvcf1="calls/selected.vcf",
        path=config['params']['annotation']['path'], #Path to the database (buildver)
        buildver=config['params']['annotation']['buildver']
    output:
        o1="annotated/annotated_all_snps_"+config['params']['annotation']['output']
    resources:
        mem_mb=config['mem_mb']
    threads: config['threads']
    conda:
        "envs/perl.yaml"
    shell:
        "perl ./rules/scripts/table_annovar.pl {input.gvcf1} {input.path} -buildver {input.buildver} -out {input.i3} -remove -protocol avsnp151 -operation g -nastring . -vcfinput"
