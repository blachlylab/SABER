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

def umi_switch(wildcards):
    if(config["use_umi"]):
        return "output/fastq_umi/"+wildcards.sample+".fastq"
    else:
        return "output/cutadapt3p/"+wildcards.sample+".fastq"

rule all:
    input:"analysis.done"

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
    output:"output/fastqc/{sample}.html"
    conda:"envs/tools-env.yaml"
    shell:"fastqc -o fastqc/ {input}"

#rule multiqc:
    #input:
# run multiqc to summarize the qc files
#system(paste0('multiqc -d ',output_folder," -o ",output_folder))

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

rule altersam:
    input:
        sam="output/needleall/{sample}.sam",
        header="references/sam_header.tsv"
    output:"output/bam/{sample}.bam"
    conda:"envs/tools-env.yaml"
    shell:"cat {input.header} <(grep -v ^@ {input.sam} | awk '{{if ($4==1){{print}}}}') | samtools view -hb | samtools sort> {output}"

rule index:
    input:"output/bam/{sample}.bam"
    output:"output/bam/{sample}.bam.bai"
    conda:"envs/tools-env.yaml"
    shell:"samtools index {input}"

rule analysis:
    input:expand("output/bam/{sample}.bam.bai",sample=SAMPLES)
    output:touch("analysis.done")
    conda:"envs/r-env.yaml"
    shell:"""
    Rscript scripts/installrpy.R
    Rscript scripts/analysis.R
    """
    

# fix sam file header and select only reads matching at position 1
# for (i in 1:length(sam_temp_files)) {
#   samdf<-read.delim(sam_temp_fp[i], sep="\t", header = FALSE)#read in sam file
#   samdf<-samdf[-(1:2),]#chop off the old header
#   samdf<-samdf[which(samdf$V4=="1"),]#select only reads mapping to coordinate 1
#   sam_header<-read.table("references/sam_header.csv", fill=TRUE, header=FALSE, sep=",", colClasses=(rep("character",13)))# read in standard sam header
#   names(sam_header)<-paste("V", 1:13, sep="")
#   samdf<-rbind(sam_header,samdf)
#   write_tsv(samdf, na = "", path = sam_temp_fp[i],col_names = FALSE, append=FALSE)
# }
# sam_temp_files<-list.files(sam_temp); sam_temp_fp<-paste0(sam_temp,"/",sam_temp_files)

# convert sam to bam
# foreach(i=1:length(sam_temp_files)) %dopar% {
#   cmd<-paste0("samtools view -S -b ",sam_temp_fp[i]," > ",bam_temp,"/",mid_names[i],".bam")
#   message(cmd,"\n"); system(cmd)
# }

# bam_temp_files<-list.files(path = bam_temp); bam_temp_fp<-paste0(bam_temp,"/",bam_temp_files)


# sort and index bam
# foreach(i=1:length(bam_temp_fp)) %dopar% {
#   cmd<-paste0("samtools sort ",bam_temp_fp[i]," -o ",bam_temp_fp[i])
#   message(cmd, "\n"); system(cmd)
# }

# foreach(i=1:length(bam_temp_fp)) %dopar% {
#   cmd<-paste0("samtools index ",bam_temp_fp[i])
#   message(cmd, "\n"); system(cmd)
# }

# system(paste0("cp -r temp/bam_temp ",bam_output_folder))#move final bam files and indices to output

#deduplicate, sort and index UMI tagged reads (optional)
# if (use_UMI==TRUE){
#   foreach(i=1:length(fastq_umi_out_files)) %dopar% {
#     cmd<-paste0("umi_tools dedup --method=unique -I ",bam_temp_fp[i]," --output-stats=",bam_dedup_output_folder,"/",mid_names[i]," -S ",bam_dedup,"/",mid_names[i],".dedup.bam") 
#     message(cmd, "\n"); system(cmd)
#   }
#   bam_dedup_files<-list.files(path = bam_dedup); bam_dedup_fp<-paste0(bam_dedup,"/",bam_dedup_files)
  
#   foreach(i=1:length(bam_dedup_fp)) %dopar% {
#     cmd<-paste0("samtools sort ",bam_dedup_fp[i]," -o ",bam_dedup_fp[i])
#     message(cmd, "\n"); system(cmd)
#   }
  
#   foreach(i=1:length(bam_dedup_fp)) %dopar% {
#     cmd<-paste0("samtools index ",bam_dedup_fp[i])
#     message(cmd, "\n"); system(cmd)
#   }
#   system(paste0("cp -r temp/bam_dedup ", bam_dedup_output_folder))
# }

#build metadata for experiment
# if (use_UMI==TRUE){
#   bam_fnames<-bam_dedup_fp
# } else {
#   bam_fnames <- bam_temp_fp
# }
# group_desig<-rep(group_name, times=length(bam_fnames))
# md<-read.csv("references/blank_metadata.csv", header = TRUE)
# newrow<-data.frame(bamfile=bam_fnames, directory=getwd(),Short.name=mid_names,Targeting.type="",sgRNA1="",sgRNA2="",Group=group_desig)
# md<-rbind(md,newrow)
