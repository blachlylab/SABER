import glob
import os.path
import re
configfile: "config.json"
#load samples from fastq folder
#samples<-list.files(paste0(getwd(),"/fastq")) # loads all files in the fastq folder
#samples_fp<-paste0(getwd(),"/fastq/",samples)
fq_files=glob.glob(os.path.join(config['fastq_dir'], "*.fastq.gz"))
r = re.compile("(.*)(_S[0-9].*)?(_L00[1-9])?_R[12](_001)?.fastq.gz")
fq_files= list(filter(r.match,fq_files))
print(fq_files)
SAMPLES=list(set([r.match(x).group(1)[len(config['fastq_dir'])+1::] for x in fq_files]))
#SAMPLES, = glob_wildcards(os.path.join(config['fastq_dir'], "{sid}(_S[0-9].*)?(_L00[1-9])?_R[12](_001)?.fastq.gz"))
print(SAMPLES)

def remap(wildcards):
    # Miseq files end in _001.fastq.gz ; NWCH Hiseq files do not
    # r1list = glob.glob(os.path.join(config['fastq_dir'], wildcards.sample + "_*_R1.fastq.gz")) # usu ends in _001.fastq.gz; not NWCH
    # r2list = glob.glob(os.path.join(config['fastq_dir'], wildcards.sample + "_*_R2.fastq.gz"))
    r1_glob_pattern = os.path.join(config['fastq_dir'], wildcards.sample + "_*R1")
    r2_glob_pattern = os.path.join(config['fastq_dir'], wildcards.sample + "_*R2")
    r1_glob_pattern += "*.fastq.gz"
    r2_glob_pattern += "*.fastq.gz"

    r1list = glob.glob(r1_glob_pattern)
    r2list = glob.glob(r2_glob_pattern)

    r1file = r1list[0]
    r2file = r2list[0]
    return [r1file, r2file]

def umi_trim_switch(wildcards):
    if(config["use_umi"]):
        return "output/fastq_umi/"+wildcards.sample+".fastq"
    else:
        return "output/cutadapt3p/"+wildcards.sample+".fastq"

def umi_bam_switch(wildcards):
    if(config["use_umi"]):
        return expand("output/dedup/{sample}.bam.bai",sample=SAMPLES)
    else:
        return expand("output/bam/{sample}.bam.bai",sample=SAMPLES)


rule all:
    input:"analysis.done","output/multiqc/multiqc_report.html"

"""
NOT USED
"""
rule decompress:
    input:"fastq/{sample}.fastq.gz"
    output:"fastq/{sample}.fastq"
    # TODO: should I assume gzip is installed?
    shell:"gzip -d -c {input} > {output}"

rule merge_read_pairs:
    input:remap
    output:"output/fastq_merged/{sample}.assembled.fastq"
    log:"logs/pear/{sample}.PEARreport.txt"
    conda:"envs/tools-env.yaml"
    threads:1
    shell:"pear -j {threads} -f {input[0]} -r {input[1]} -o output/fastq_merged/{wildcards.sample} > {log}"

rule trimmomatic:
    input:"output/fastq_merged/{sample}.assembled.fastq"
    output:"output/fastq_trimmed/{sample}.trimmed.fastq"
    conda:"envs/tools-env.yaml"
    shell:"trimmomatic SE {input} {output} SLIDINGWINDOW:4:15 MINLEN:100"

#run cutadapt to keep only fastqs that have the full flanking primer sequences
#5 prime adapter = V6/7F:  TCGAGCTCAAGCTTCGG
rule cutadapt5p:
    input:"output/fastq_trimmed/{sample}.trimmed.fastq"
    output:"output/cutadapt5p/{sample}.fastq"
    conda:"envs/tools-env.yaml"
    #changed to allow full adapter sequence anywhere to accomodate UMI.  If this doesn't work change X to ^.
    shell:"cutadapt -g XTCGAGCTCAAGCTTCGG --discard-untrimmed -e 0.01 --action=none -o {output} {input}"

#3 prime adapter = V6/7R:  GACCTCGAGACAAATGGCAG (reverse complement of the primer sequence 5'-3')    
rule cutadapt3p:
    input:"output/cutadapt5p/{sample}.fastq"
    output:"output/cutadapt3p/{sample}.fastq"
    conda:"envs/tools-env.yaml"
    #changed to allow full adapter sequence anywhere to accomodate UMI.  If this doesn't work change X to ^.
    shell:"cutadapt -a GACCTCGAGACAAATGGCAG$ --discard-untrimmed -e 0.01 --action=none -o {output} {input}"

# run fastqc on trimmed and demultiplexed input files
rule fastqc:
    input:"output/cutadapt3p/{sample}.fastq"
    output:"output/fastqc/{sample}_fastqc.html","output/fastqc/{sample}_fastqc.zip"
    conda:"envs/tools-env.yaml"
    shell:"fastqc -o output/fastqc/ {input}"

# run multiqc to summarize the qc files
rule multiqc:
    input:expand("output/fastqc/{sample}_fastqc.html",sample=SAMPLES)
    output:"output/multiqc/multiqc_report.html"
    conda:"envs/tools-env.yaml"
    shell:"multiqc -d output/fastqc/ -o output/multiqc/"

# extract UMI tags (optional)
rule umi_extract:
    input:"output/cutadapt3p/{sample}.fastq"
    output:"output/fastq_umi/{sample}.fastq"
    log:"logs/umi_extract/{sample}.log"
    conda:"envs/tools-env.yaml"
    shell:"umi_tools extract --stdin={input} --bc-pattern=NNNNNNNNNN --log={log} --stdout={output}"

rule needleall:
    input:umi_switch
    output:"output/needleall/{sample}.sam"
    log:"logs/needleall/{sample}.error"
    conda:"envs/tools-env.yaml"
    shell:"needleall -aformat3 sam -gapextend 0.25 -gapopen 10.0 -awidth3=5000 -asequence references/gestalt_pipeline4.fasta -bsequence {input} -outfile {output} -errfile {log}"

# fix sam file header and select only reads matching at position 1
rule altersam:
    input:
        sam="output/needleall/{sample}.sam",
        header="references/sam_header.tsv"
    output:"output/bam/{sample}.bam"
    conda:"envs/tools-env.yaml"
    shell:"cat {input.header} <(grep -v ^@ {input.sam} | awk '{{if ($4==1){{print}}}}') | samtools view -hb | samtools sort> {output}"

rule index:
    input:"{fn}.bam"
    output:"{fn}.bam.bai"
    conda:"envs/tools-env.yaml"
    shell:"samtools index {input}"

rule dedup:
    input:
        bam="output/bam/{sample}.bam",
        bai="output/bam/{sample}.bam.bai"
    output:"output/dedup/{sample}.bam"
    conda:"envs/tools-env.yaml"
    shell:"umi_tools dedup --method=unique -I {input.bam} --output-stats= output/dedup/{wildcards.sample} -S {output}"


rule analysis:
    input:umi_bam_switch
    output:touch("analysis.done")
    conda:"envs/r-env.yaml"
    shell:"""
    Rscript scripts/installrpy.R
    Rscript scripts/analysis.R
    """
    
