import os
import pandas as pd
import fnmatch

samples = []
with open("../sra_number", "r") as file:
    for line in file:
        if line != ''
        samples.append(line.strip())
print(samples)


rule all:
    input: 
        expand("/home/zzy/graduation_project/data/{sample}.gene.count.txt",sample=samples),
        expand("/home/zzy/graduation_project/data/{sample}.gene.count.txt.summary",sample=samples)
    output:
        expand("/home/zzy/graduation_project/data/{sample}.rawcount.tsv", sample=samples)
    run:
        exp = pd.DataFrame()
        for i in input:
            tab = pd.read_csv(i, skiprows=1, sep='\t')
            exp[ i.split('/')[-1].split('.')[0] ] = tab.iloc[:, 6]
        else:
            exp.index = tab['Geneid']
        exp.to_csv(output[0], sep='\t')
        
rule fastq_dump:
    input: expand("{sample}.sra",sample=samples)
    output: 
        "/home/zzy/graduation_project/data/{sample}_1.fastq",
        "/home/zzy/graduation_project/data/{sample}_2.fastq"
    threads: 8
    shell: 
        """
        cd /home/zzy/graduation_project/data/tmp
        /home/zzy/graduation_project/sratoolkit/bin/fastq-dump --split-3 {input} 
        mv {sample}_1.fastq {output[0]}
        mv {sample}_2.fastq {output[1]}
        """

rule fastp:
    input: rules.fastq_dump.output
    output: 
        "/home/zzy/graduation_project/data/{sample}_R1.clean.fq",
        "/home/zzy/graduation_project/data/{sample}_R2.clean.fq",
    threads: 8
    shell:
        """
        fastp -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} --detect_adapter_for_pe \
        -q 30 -e 30 -L 30 -c -g \
        -M 1 --cut_fount 1 --cut_tail 1\
        -w {threads}
        """

rule STAR:
    input: rules.fastp.output 
    output: 
        "/home/zzy/graduation_project/data/{sample}.Aligned.sortedByCoord.out.bam",
        temp("/home/zzy/graduation_project/data/{sample}.ReadsPerGene.out.tab"),
        "/home/zzy/graduation_project/data/{sample}.Log.final.out",
        temp("/home/zzy/graduation_project/data/{sample}.Log.out"),
        temp("/home/zzy/graduation_project/data/{sample}.SJ.out"),
        temp("/home/zzy/graduation_project/data/{sample}.Log.progress.out"),
    threads: 16
    params:
        STAR="/home/chensy/anaconda3/envs/snakemake/bin/STAR",
        genome="/workspace/chensy/Immunity/0.Script/reference/GRCh38_STAR_index"
    shell:
        """
        {params.STAR} --genomeDir {params.genome} --readFilesIn {input} --runThreadN {threads} \
        --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 --alighIntronMin 20 --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 ==outSAMunmapped Within --outFilterType BySJout \
        --outFilterNamePrefix /home/zzy/graduation_project/data/{wildcards.sample} \
        --outSAMstrandField intronMotif \
        --quantMode GeneCounts \
        --outSAMtype BAM SortedByCoordinate --sjdbScore 1
        """ 

rule featureCounts:
    input: "/home/zzy/graduation_project/data/{sample}.Aligned.sortedByCoord.out.bam"
    output: 
        "/home/zzy/graduation_project/data/{sample}.gene.count.txt",
        "/home/zzy/graduation_project/data/{sample}.gene.count.txt.summary"
    threads: 16
    shell:
        "featureCounts -T {threads} -t exon -g gene_id -a /workspace/chensy/Immunity/0.Script/reference/hg38.refGene.gtf -p -o {output[0]} {input}"


