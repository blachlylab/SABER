import glob
import os.path
import re
configfile: "config.json"
#load samples from fastq folder
#samples<-list.files(paste0(getwd(),"/fastq")) # loads all files in the fastq folder
#samples_fp<-paste0(getwd(),"/fastq/",samples)
# fq_files=glob.glob(os.path.join(config['fastq_dir'], "*.fastq.gz"))
# r = re.compile("(.*)(_S[0-9].*)?(_L00[1-9])?_R[12](_001)?.fastq.gz")
# fq_files= list(filter(r.match,fq_files))
# print(fq_files)
# SAMPLES=list(set([r.match(x).group(1)[len(config['fastq_dir'])+1::] for x in fq_files]))
# #SAMPLES, = glob_wildcards(os.path.join(config['fastq_dir'], "{sid}(_S[0-9].*)?(_L00[1-9])?_R[12](_001)?.fastq.gz"))
# print(SAMPLES)

def INPUTFOLDER(wildcards):
    output="fastq"
    if wildcards.project == "testing":
        output="testing_"+output
    elif wildcards.project == "validation":
        output="validation_"+output
    return output

def INPUTSAMPLES(wildcards):
    output = INPUTFOLDER(wildcards)
    fq_files=glob.glob(os.path.join(output, "*.fastq.gz"))
    r = re.compile("(.*)(_S[0-9].*)?(_L00[1-9])?_R[12](_001)?.fastq.gz")
    fq_files= list(filter(r.match,fq_files))
    SAMPLES=list(set([r.match(x).group(1)[len(output)+1::] for x in fq_files]))
    print(SAMPLES)
    return SAMPLES

def remapSampleToFastq(wildcards):
    # Miseq files end in _001.fastq.gz ; NWCH Hiseq files do not
    # r1list = glob.glob(os.path.join(config['fastq_dir'], wildcards.sample + "_*_R1.fastq.gz")) # usu ends in _001.fastq.gz; not NWCH
    # r2list = glob.glob(os.path.join(config['fastq_dir'], wildcards.sample + "_*_R2.fastq.gz"))
    r1_glob_pattern = os.path.join(INPUTFOLDER(wildcards), wildcards.sample + "_*R1")
    r2_glob_pattern = os.path.join(INPUTFOLDER(wildcards), wildcards.sample + "_*R2")
    r1_glob_pattern += "*.fastq.gz"
    r2_glob_pattern += "*.fastq.gz"

    r1list = glob.glob(r1_glob_pattern)
    r2list = glob.glob(r2_glob_pattern)

    r1file = r1list[0]
    r2file = r2list[0]
    return [r1file, r2file]

def umi_trim_switch(wildcards):
    if(config["use_umi"]):
        return wildcards.project+"_output/fastq_umi/"+wildcards.sample+".fastq"
    else:
        return wildcards.project+"_output/cutadapt3p/"+wildcards.sample+".fastq"

def umi_bam_switch(wildcards):
    if(config["use_umi"]):
        return expand(wildcards.project+"_output/dedup/{sample}.bam.bai",sample=INPUTSAMPLES(wildcards))
    else:
        return expand(wildcards.project+"_output/bam/{sample}.bam.bai",sample=INPUTSAMPLES(wildcards))

def multiqc_input(wildcards):
    return expand("{project}_output/fastqc/{sample}_fastqc.html",sample=INPUTSAMPLES(wildcards),project=wildcards.project)


rule normal_run:
    input: "gestalt_output/analysis.done","gestalt_output/multiqc/multiqc_report.html"

rule testing:
    input: "testing_output/analysis.done","testing_output/multiqc/multiqc_report.html"

rule validation:
    input: "validation_output/analysis.done","validation_output/multiqc/multiqc_report.html"


"""
NOT USED
"""
rule decompress:
    input:"fastq/{sample}.fastq.gz"
    output:"fastq/{sample}.fastq"
    # TODO: should I assume gzip is installed?
    shell:"gzip -d -c {input} > {output}"

rule merge_read_pairs:
    input:remapSampleToFastq
    output:"{project}_output/fastq_merged/{sample}.assembled.fastq"
    log:"{project}_logs/pear/{sample}.PEARreport.txt"
    conda:"envs/tools-env.yaml"
    threads:1
    shell:"pear -j {threads} -f {input[0]} -r {input[1]} -o output/fastq_merged/{wildcards.sample} > {log}"

rule trimmomatic:
    input:"{project}_output/fastq_merged/{sample}.assembled.fastq"
    output:"{project}_output/fastq_trimmed/{sample}.trimmed.fastq"
    conda:"envs/tools-env.yaml"
    shell:"trimmomatic SE {input} {output} SLIDINGWINDOW:4:15 MINLEN:100"

#run cutadapt to keep only fastqs that have the full flanking primer sequences
#5 prime adapter = V6/7F:  TCGAGCTCAAGCTTCGG
rule cutadapt5p:
    input:"{project}_output/fastq_trimmed/{sample}.trimmed.fastq"
    output:"{project}_output/cutadapt5p/{sample}.fastq"
    params:
        adapter= "XCGCAGAGAGGCTCCGTG" if config["use_umi"] else "XTCGAGCTCAAGCTTCGG" 
    conda:"envs/tools-env.yaml"
    #changed to allow full adapter sequence anywhere to accomodate UMI.  If this doesn't work change X to ^.
    shell:"cutadapt -g {params.adapter} --discard-untrimmed -e 0.01 --action=none -o {output} {input}"

#3 prime adapter = V6/7R:  GACCTCGAGACAAATGGCAG (reverse complement of the primer sequence 5'-3')    
rule cutadapt3p:
    input:"{project}_output/cutadapt5p/{sample}.fastq"
    output:"{project}_output/cutadapt3p/{sample}.fastq"
    conda:"envs/tools-env.yaml"
    #changed to allow full adapter sequence anywhere to accomodate UMI.  If this doesn't work change X to ^.
    shell:"cutadapt -a GACCTCGAGACAAATGGCAG$ --discard-untrimmed -e 0.01 --action=none -o {output} {input}"

# run fastqc on trimmed and demultiplexed input files
rule fastqc:
    input:"{project}_output/cutadapt3p/{sample}.fastq"
    output:"{project}_output/fastqc/{sample}_fastqc.html","{project}_output/fastqc/{sample}_fastqc.zip"
    conda:"envs/tools-env.yaml"
    shell:"fastqc -o output/fastqc/ {input}"

# run multiqc to summarize the qc files
rule multiqc:
    input: multiqc_input
    output:"{project}_output/multiqc/multiqc_report.html"
    conda:"envs/tools-env.yaml"
    shell:"multiqc -d {wildcards.project}_output/fastqc/ -o {wildcards.project}_output/multiqc/"

# extract UMI tags (optional)
rule umi_extract:
    input:"{project}_output/cutadapt3p/{sample}.fastq"
    output:"{project}_output/fastq_umi/{sample}.fastq"
    log:"{project}_logs/umi_extract/{sample}.log"
    conda:"envs/tools-env.yaml"
    shell:"umi_tools extract --stdin={input} --bc-pattern=XXXXXXXXXXXXXXXXXNNNNNNNNNN --log={log} --stdout={output}"

rule needleall:
    input:umi_trim_switch
    output:"{project}_output/needleall/{sample}.sam"
    log:"{project}_logs/needleall/{sample}.error"
    conda:"envs/tools-env.yaml"
    shell:"needleall -aformat3 sam -gapextend 0.25 -gapopen 10.0 -awidth3=5000 -asequence references/gestalt_pipeline4.fasta -bsequence {input} -outfile {output} -errfile {log}"

# fix sam file header and select only reads matching at position 1
rule altersam:
    input:
        sam="{project}_output/needleall/{sample}.sam",
        header="references/sam_header.tsv"
    output:"{project}_output/bam/{sample}.bam"
    conda:"envs/tools-env.yaml"
    shell:"cat {input.header} <(grep -v ^@ {input.sam} | awk '{{if ($4==1){{print}}}}') | samtools view -hb | samtools sort> {output}"

rule index:
    input:"{fn}.bam"
    output:"{fn}.bam.bai"
    conda:"envs/tools-env.yaml"
    shell:"samtools index {input}"

rule dedup:
    input:
        bam="{project}_output/bam/{sample}.bam",
        bai="{project}_output/bam/{sample}.bam.bai"
    output:"{project}_output/dedup/{sample}.bam"
    conda:"envs/tools-env.yaml"
    shell:"umi_tools dedup --method=unique -I {input.bam} --output-stats= output/dedup/{wildcards.sample} -S {output}"


rule analysis:
    input:umi_bam_switch
    output:touch("{project}_output/analysis.done")
    conda:"envs/r-env.yaml"
    params: 
        script= lambda wildcards: "analysis" if wildcards.project=="testing" else "validation",
        bams= lambda wildcards: wildcards.project+"_output/dedup" if config["use_umi"] else wildcards.project+"_output/bam"
    shell:"""
    Rscript scripts/installrpy.R
    Rscript scripts/{params.script}.R {params.bams}
    """

rule testing_analysis:
    input:umi_bam_switch
    output:touch("testing_output/analysis.done")
    conda:"envs/r-env.yaml"
    shell:"""
    Rscript scripts/installrpy.R
    Rscript scripts/analysis.R
    """

rule analysis_validation:
    input:umi_bam_switch
    output:touch("validation_output/validation.done")
    conda:"envs/r-env.yaml"
    shell:"""
    Rscript scripts/installrpy.R
    Rscript scripts/validation.R
    """

