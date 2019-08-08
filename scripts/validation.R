args=commandArgs(trailingOnly=TRUE)
####libraries####
library(CrispRVariants) 
library("Rsamtools")
library(rtracklayer)
library(GenomicFeatures)
library(gridExtra)
library(GenomicRanges)
library(RColorBrewer)
library(rPython)
library(sqldf)
library(rmarkdown)
library(foreach)
library(doParallel)
library(stringr)
library(cowplot)
library(reshape2)
library(RColorBrewer)
library(ggdendro)
library(cowplot)
library(scales)
library(viridis)
library(ggExtra)
library(tidyverse)
library(ggrepel)


####version notes####
# version 3.0 uses paired illumina reads from genewiz and needleall aligner\
# version 3.1 changes thresholding method for stacked bar plots and csvs to proportion of all reads
# version 3.2 parallelizes needleall
# no umis
# You need a genome file in a folder named "genome" and your fastqs in a folder named "fastq" for this version to work
# version 4.0 combines main pipeline with clone comparer
# version 4.1 spits out fraction informative data too
# version 4.2_normal_samples is designed to eval normal samples
# analysis_v1.0 derived from pipeline v4.2
# validation_v1.0 uses optimized settings from analysis v1.0 and adds in side by side comparison between groups
# validation_v2.0 cleans up v1.0

####general functions####
se <- function(x) sqrt(var(x, na.rm=TRUE)/length(x))

data_summary<-function(x){
  m<-median(x)
  ymin<-m-se(x)
  ymax<-m+se(x)
  return(c(y=m, ymin=ymin,ymax=ymax))
}

data_summary1<-function(x){
  m<-median(x)
  ymin<-m-(IQR(x)/2)
  ymax<-m+(IQR(x)/2)
  return(c(y=m, ymin=ymin,ymax=ymax))
}

data_summary2<-function(x){
  m<-mean(x)
  ymin<-m-se(x)
  ymax<-m+se(x)
  return(c(y=m, ymin=ymin,ymax=ymax))
}

data_summary3<-function(x){
  m<-mean(x)
  ymin<-m-sd(x)
  ymax<-m+sd(x)
  return(c(y=m, ymin=ymin,ymax=ymax))
}

# diversity index functions
shannon_func<-function(x){
  prop<-x/sum(x)
  shannon_entropy<--1*sum(prop*log2(prop))
  return(shannon_entropy)
}

simpson_func<-function(x){
  prop<-x/sum(x)
  simpson_diversity<-1/sum(prop^2)
  return(simpson_diversity)
}

general_diversity_func<-function(x,q){
  prop<-x/sum(x)
  general_diversity<-(sum(prop^q))^(1/(1-q))
  return(general_diversity)
}

####input data####
#Be sure this is correct before running!!!                                                
genome<-"gestalt_pipeline4.fasta" # include .fasta.  Genome file has to be in genome folder.
genome_fp<-paste0(getwd(),"/references/",genome)
output_file<-"training_vs_validation_set"#for the cripsrvariants plot and prefix on sample output csv and vaf plot
group_name<-"training_vs_validation_set"

threshold<-200# number of reads below which they don't appear on the big crispr plot.  THis only affects the big cripsrplot!

thresh<-0;tvaf<-1;tpaf<-0.1#threshold settings for main analysis pipeline

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

# make working subdirectories and create variables 
output_folder<-"validation_output"
bam_files<-list.files(path = args[2],".*.bam")
bam_fnames<-paste0(args[2],"/",bam_files)
mid_names<-substr(bam_files, 6,10)

group_desig<-rep(group_name, times=length(bam_fnames))
md<-read.csv("references/blank_metadata.csv", header = TRUE)
newrow<-data.frame(bamfile=bam_fnames, directory=getwd(),Short.name=mid_names,Targeting.type="",sgRNA1="",sgRNA2="",Group=group_desig)
md<-rbind(md,newrow)

#create target region
gd <- rtracklayer::import("references/gestalt2.bed")
gdl <- GenomicRanges::resize(gd, width(gd) + 0, fix = "center") #resize region for analysis
reference0<-read_file("references/gestalt2_ref.txt")
reference1<-substr(reference0,1,310)
reference<-Biostrings::DNAString(reference1)
reference

# make the crispr set
crispr_set <- readsToTarget(bam_fnames, 
                            target = gd, 
                            reference = reference, 
                            names = md$Short.name, 
                            target.loc = 16, 
                            collapse.pairs = FALSE,
                            split.snv=FALSE)#split.snv=FALSE adds SNVs into no variant count
# plot the variants
ps<-37
pam_seq<-seq(ps,280,27)

while (!is.null(dev.list())) dev.off()
p <- plotVariants(crispr_set, 
                  col.wdth.ratio = c(1,1),
                  plotAlignments.args = list(pam.start = pam_seq, #c(37,64), #draws a 3-nt box starting including the position noted
                                             target.loc = pam_seq-3, #draws a vertical line after the position noted
                                             guide.loc = IRanges::IRanges(pam_seq-20,pam_seq+2), #first parameter - beginning of target sequence, second - end of target sequence
                                             min.count = threshold,
                                             tile.height = 0.9),
                  plotFreqHeatmap.args = list(min.count = threshold,
                                              plot.text.size = 3, 
                                              x.size = 8, 
                                              group = group_desig,
                                              legend.text.size = 8,
                                              legend.key.height = grid::unit(0.5, "lines")))
dev.copy2pdf(file=paste0(output_folder,"/",output_file,".pdf"), width = 36, height = 36)  #for 10 samples use 24 x 24, for 30-40 samples use 48 x 48

#####generate informative, no variant and common variant tables with vaf and paf thresholds####
generate_v1<-function(x,thresh,vaf_thresh,paf_thresh,vc_prop,output_folder,output_file){
  df<-x#isolate df 
  rnames<-row.names(df)#duplicate row names as new variable
  df$rnames<-rnames#add row name variabe to df
  df<-df[which(df$rnames!="Other"),]# removes other
  common_vars<-names(which(rowSums(vc_prop>vaf_thresh)/ncol(vc_prop)>(paf_thresh)))#returns the names of rows with a propotion greater than vaf_thresh occuring in more than paf_thresh
  remove<-c("no variant")# removes only no variant.  Leaves "Other" in.  Maybe need to pull it out.
  common_vars<-setdiff(common_vars,remove)
  for (i in 1:length(common_vars)) df$rnames<-gsub(paste0("^",common_vars[i],"$"),"common variant",df$rnames)#replaces all common variants as "common variant"
  #df<-df[which(df[1]/sum(df[1])>thresh | df$rnames=="no variant"),]# removes reads below threshold and keeps no variant##20190509 changed to use proportional threshold
  common_variant_sum<-sum(df[1][which(df$rnames=="common variant"),])#sums all common variant calls
  df$is_common_variant<-NA#adds a column that will hold the common variant marker
  df[nrow(df) + 1,] = list(common_variant_sum,"common variant sum","yes")#adds a new row to df with summed novariant calls.  Yes is a no variant marker
  df$is_no_variant<-NA#adds a column that will hold the no_variant marker
  df[which(df$rnames == "no variant"),4]<-"yes"
  df<-df[which(df$rnames!="common variant"),,drop=FALSE]#drop all constituents of no variant sum
  df<-df[which(df$rnames!="Other"),,drop=FALSE]#drop other reads
  df<-df[order(df$is_no_variant, df$is_common_variant,-df[1]),]# sort on no variant marker then on read counts
  df<-df[which(df[1]/sum(df[1])>thresh | df$is_no_variant=="yes" | df$is_common_variant=="yes"),]# removes reads below read count threshold, keeping no variant calls and common variant calls
  n_variants<-length(df$rnames)#total number of reads in each specimen
  total_reads<-sum(df[1])#sum of reads for each specimen
  allele_freq<-df[1]/total_reads# calculate allele freq by line
  names(allele_freq)[1] <- "allele_freq"# rename column for allele_freq
  df<-cbind(df,allele_freq)#add allele_freq column to the working df
  specimen_id<-names(df)[1]# pull the specimen name from the column title
  specimen<-rep_len(specimen_id, n_variants)# make a new column for the specimen id
  df<-cbind(df,specimen)# bind it to the working df
  #fix levels
  index<-as.factor(c(1:length(df$rnames)))#create index string as long as the number of variants, make it as a factor
  index<-factor(index, levels=rev(levels(index)))#reverse the levels
  df<-cbind(df,index)#bind it to the working df
  dir.create(paste0(output_folder,"/thresh_",thresh,"_tvaf_",vaf_thresh,"_tpaf_",paf_thresh))
  write.csv(df,file = paste0(output_folder,"/thresh_",thresh,"_tvaf_",vaf_thresh,"_tpaf_",paf_thresh,"/",output_file,"_",specimen_id,".csv"))
}

vc_all <- as.data.frame(variantCounts(crispr_set), stringsAsFactors=FALSE)#big data frame of read counts
vc_list<-list()
for(i in 1:length(mid_names)){vc_list[[length(vc_list)+1]]<-vc_all[i]}#puts individual samples into a list of vectors
vc_prop<-variantCounts(crispr_set, result = "proportions")# generates a matrix of variant allele fractions, similar to vc_all

parLapply(cl,vc_list,generate_v1,thresh = thresh, vaf_thresh = tvaf, paf_thresh = tpaf,vc_prop = vc_prop,output_folder = output_folder, output_file = output_file)

####identify informative and uninformative samples using kmeans clustering####
#load the data
infiles<-list.files(paste0(output_folder,"/thresh_",thresh,"_tvaf_",tvaf,"_tpaf_",tpaf), pattern = "\\.csv", full.names = TRUE)

fidf_0<-read.csv(infiles[1])#reads in variant csvs for each sample
fi<-1-sum(c(fidf_0[1,2],fidf_0[2,2]))/sum(fidf_0[2])#calculates fraction informative from read counts
fidf_1<-fidf_0[,c(7,2,3)]#keeps only the columns we want
colnames(fidf_1)<-c("specimen", "counts", "variant")#renames the columns
fidf_1$FI<-rep(x = fi,times = nrow(fidf_0))#makes a new column holding fraction informative
for (i in 2:length(infiles)){#for loop to read in all of the samples analyzed
  newdf<-read.csv(infiles[i])
  fi<-1-sum(c(newdf[1,2],newdf[2,2]))/sum(newdf[2])
  newdf<-newdf[,c(7,2,3)]
  colnames(newdf)<-c("specimen", "counts", "variant")
  newdf$FI<-rep(x = fi,times = nrow(newdf))
  fidf_1<-rbind(fidf_1,newdf)
}
fidf_2<-fidf_1[,c(1,4)]#gets rid of the variant-level data for each sample since you don't need it here
fidf_3<-unique(fidf_2)#gets rid of duplicate rows
fidf_3$rank<-rank(-(fidf_3$FI), ties.method = "first")#ranks samples based on fraction informative with 1 = most informative

#kmeans
km<-kmeans(fidf_3$FI, centers = 2)#generates kmeans of 2 clusters based on fraction informative
fidf_3$cluster<-as.factor(km$cluster)#puts cluster vector into the fi dataframe as a new column


top_clust<-fidf_3[which(fidf_3$rank==1),4]#identifies the top cluster since kmeans function assigns cluster 1 at random
bottom_clust<-fidf_3[which(fidf_3$rank==length(fidf_3$rank)),4]#identifies bottom cluster
fidf_3$cluster<-as.character(fidf_3$cluster)
fidf_3$cluster[fidf_3$cluster==top_clust]<-"Informative"#relabel clusters so they are consistent
fidf_3$cluster[fidf_3$cluster==bottom_clust]<-"Not Informative"#ditto


md<-read.csv(list.files("references", pattern = "clone_md.csv", full.names = TRUE))
colnames(md)<-c("specimen", "class", "batch","t_v")
fidf_4<-merge(fidf_3,md,by = "specimen", all.x = TRUE)


#plot the fi distribution with kmeans
fi_dist<-ggplot(fidf_3, aes(x=factor(fidf_3$rank), y=fidf_3$FI, fill=fidf_3$cluster))
fi_dist<-fi_dist+
  geom_bar(stat = "identity", color = "black")+
  scale_fill_manual(values = c("#72ce55","#fde725"),name = "")+
  scale_x_discrete(breaks = c(5,10,15,20,25,30,35,40,45,50,55,60))+
  theme_cowplot(font_size = 11)+
  theme(legend.position = c(0.6,0.95))+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8))+
  labs(x="Sample Index", y="Fraction Informative")
fi_dist
save_plot(plot = fi_dist, filename = "output_validation/fidist_f23a1.pdf", base_height = 3, base_width = 4)

#batch plot
fi_batch<-ggplot(fidf_4, aes(x=factor(fidf_4$rank), y=fidf_4$FI, fill=fidf_4$t_v))
fi_batch<-fi_batch+
  geom_bar(stat = "identity", color = "black")+
  scale_x_discrete(breaks = c(5,10,15,20,25,30,35,40,45,50,55,60))+
  scale_fill_manual(values = c("#DC0000","#3C5488"),
                    name = "",
                    breaks=c("training","validation"),
                    labels=c("Training","Validation"))+
  theme_cowplot(font_size = 11)+
  theme(legend.position = c(0.6,0.95))+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8))+
  labs(x="Sample Index", y="Fraction Informative")
fi_batch
save_plot(plot = fi_batch, filename = "output_validation/fibatch_7f5f3.pdf", base_height = 3, base_width = 4)

all_samples<-fidf_4$specimen
informative_samples<-filter(fidf_4, fidf_4$cluster=="Informative")

####quantifying alleles from thresholded data in informative samples####
#reassemble the long form data frame of read counts using only informative samples
inlist<-as.data.frame(cbind(as.character(all_samples), infiles))
colnames(inlist)<-c("specimen", "path")
inlist_informative1<-left_join(informative_samples,inlist, by = "specimen")
informative_path<-as.character(inlist_informative1$path)
qdf0<-read.csv(informative_path[1])
colnames(qdf0)<-c("X","read_count","rnames","is_common_variant","is_no_variant","allele_freq","specimen","index")
for (i in 2:length(informative_path)) {
  qdf_new<-read.csv(informative_path[i])
  colnames(qdf_new)<-c("X","read_count","rnames","is_common_variant","is_no_variant","allele_freq","specimen","index")
  qdf0<-rbind(qdf0,qdf_new)
}
qdf1<-qdf0[,-1]


spec_list<-unique(qdf1$specimen)

#generate a vector of shannon entropy values for the informative samples
shannon<-numeric()
for (i in 1:length(spec_list)){
  qdf_filt<-filter(qdf1, qdf1[,6]==spec_list[i])#loads in dfs one sample at a time
  qdf_filt_inf<-qdf_filt[-c(1,2),]# removes no variant and common variant
  qdf_filt_inf$index<-qdf_filt_inf$index-2#resets index count
  shannon<-c(shannon,shannon_func(qdf_filt_inf[1]))#calculates off of read count
}
shannon

#generate a vector of simpson diversity indices for the informative samples
simpson<-numeric()
for (i in 1:length(spec_list)){
  qdf_filt<-filter(qdf1, qdf1[,6]==spec_list[i])#loads in dfs one sample at a time
  qdf_filt_inf<-qdf_filt[-c(1,2),]# removes no variant and common variant
  qdf_filt_inf$index<-qdf_filt_inf$index-2#resets index count
  simpson<-c(simpson,simpson_func(qdf_filt_inf[1]))#calculates off of read count
}
simpson

#generate a vector of the number of variants over 2%
vaf2<-numeric()
for (i in 1:length(spec_list)){
  qdf_filt<-filter(qdf1, qdf1[,6]==spec_list[i])#loads in dfs one sample at a time
  qdf_filt_inf<-qdf_filt[-c(1,2),]# removes no variant and common variant
  qdf_filt_inf$index<-qdf_filt_inf$index-2#resets index count
  qdf_filt_inf$adj_allele_freq<-qdf_filt_inf$read_count/sum(qdf_filt_inf$read_count)#recalculates allele frequencies from read counts without no variant and common variant
  qdf_2percent<-filter(qdf_filt_inf, qdf_filt_inf$adj_allele_freq>=0.02)
  vaf2<-c(vaf2,nrow(qdf_2percent))
}
vaf2

#counts barcodes over 1% raw vaf
vaf1<-numeric()
for (i in 1:length(spec_list)){
  qdf_filt<-filter(qdf1, qdf1[,6]==spec_list[i])#loads in dfs one sample at a time
  qdf_filt_inf<-qdf_filt[-c(1,2),]# removes no variant and common variant
  qdf_filt_inf$index<-qdf_filt_inf$index-2#resets index count
  qdf_filt_inf$adj_allele_freq<-qdf_filt_inf$read_count/sum(qdf_filt_inf$read_count)#recalculates allele frequencies from read counts without no variant and common variant
  qdf_1percent<-filter(qdf_filt_inf, qdf_filt_inf$adj_allele_freq>=0.01)
  vaf1<-c(vaf1,nrow(qdf_1percent))
}
vaf1

#compile and summarize data
sumdf0<-as.data.frame(cbind(levels(qdf1$specimen),vaf1,vaf2,shannon,simpson))
colnames(sumdf0)<-c("Sample", "variants_over_1%", "variants_over_2%", "Shannon_Entropy", "Simpson_Diversity")
sumdf1<-left_join(sumdf0,md,by = c("Sample" = "specimen"))


sumdf_long<-gather(sumdf1, key = "attribute", value = "value_to_plot", -c(Sample,t_v,class,batch))
sumdf_long$value_to_plot<-as.numeric(sumdf_long$value_to_plot)
vaf1_long<-filter(sumdf_long, sumdf_long$attribute=="variants_over_1%")
vaf2_long<-filter(sumdf_long, sumdf_long$attribute=="variants_over_2%")
shannon_long<-filter(sumdf_long, sumdf_long$attribute=="Shannon_Entropy")
simpson_long<-filter(sumdf_long, sumdf_long$attribute=="Simpson_Diversity")

p1<-ggplot(vaf1_long, aes(x = vaf1_long$t_v, y = vaf1_long$value_to_plot))
p1<-p1+
  geom_dotplot(binaxis="y", stackdir="center", dotsize=1.5,position=position_jitter(w=0.03),aes(fill=vaf1_long$t_v)) + 
  scale_fill_manual(values=c("#DC0000","#3C5488")) + 
  scale_color_manual(values=c("black","black"))+
  stat_summary(fun.data=data_summary2, color="black", size=0.5, width=0.3, fill=c("#DC0000","#3C5488"), alpha=0.3,geom="crossbar")+
  theme_cowplot(font_size = 11)+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())+
  #scale_y_continuous(limits = c(0,5), breaks = c(0,1,2,3,4,5))+
  scale_x_discrete(breaks = c("training","validation"), labels = c("Training","Validation"))+
  labs(y="Clones with VAF>0.01")
p1
save_plot(plot = p1, filename = "output_validation/vaf_0.01_7a6ea.pdf",base_width = 2.5, base_height = 2.3)

p2<-ggplot(vaf2_long, aes(x = vaf2_long$t_v, y = vaf2_long$value_to_plot))
p2<-p2+
  geom_dotplot(binaxis="y", stackdir="center", dotsize=1.5,position=position_jitter(w=0.03),aes(fill=vaf2_long$t_v)) + 
  scale_fill_manual(values=c("#DC0000","#3C5488")) + 
  scale_color_manual(values=c("black","black"))+
  stat_summary(fun.data=data_summary2, color="black", size=0.5, width=0.3, fill=c("#DC0000","#3C5488"), alpha=0.3,geom="crossbar")+
  theme_cowplot(font_size = 11)+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10))+
  scale_x_discrete(breaks = c("training","validation"), labels = c("Training","Validation"))+
  labs(y="Clones with VAF>0.02")
p2
save_plot(plot = p2, filename = "output_validation/vaf_0.02_0dae3.pdf",base_width = 2.5, base_height = 2.3)

p3<-ggplot(shannon_long, aes(x = shannon_long$t_v, y = shannon_long$value_to_plot))
p3<-p3+
  geom_dotplot(binaxis="y", stackdir="center", dotsize=1.5,position=position_jitter(w=0.03),aes(fill=shannon_long$t_v)) + 
  scale_fill_manual(values=c("#DC0000","#3C5488")) + 
  scale_color_manual(values=c("black","black"))+
  stat_summary(fun.data=data_summary2, color="black", size=0.5, width=0.3, fill=c("#DC0000","#3C5488"), alpha=0.3,geom="crossbar")+
  theme_cowplot(font_size = 11)+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())+
  #scale_y_continuous(limits = c(0,5), breaks = c(0,1,2,3,4,5))+
  scale_x_discrete(breaks = c("training","validation"), labels = c("Training","Validation"))+
  labs(y="Shannon Entropy")
p3
save_plot(plot = p3, filename = "output_validation/shannon_c2589.pdf",base_width = 2.5, base_height = 2.3)

p4<-ggplot(simpson_long, aes(x = simpson_long$t_v, y = simpson_long$value_to_plot))
p4<-p4+
  geom_dotplot(binaxis="y", stackdir="center", dotsize=1.5,position=position_jitter(w=0.03),aes(fill=simpson_long$t_v)) + 
  scale_fill_manual(values=c("#DC0000","#3C5488")) + 
  scale_color_manual(values=c("black","black"))+
  stat_summary(fun.data=data_summary2, color="black", size=0.5, width=0.3, fill=c("#DC0000","#3C5488"), alpha=0.3,geom="crossbar")+
  theme_cowplot(font_size = 11)+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())+
  #scale_y_continuous(limits = c(0,5), breaks = c(0,1,2,3,4,5))+
  scale_x_discrete(breaks = c("training","validation"), labels = c("Training","Validation"))+
  labs(y="Simspon Diversity")
p4
save_plot(plot = p4, filename = "output_validation/simpson_a6476.pdf",base_width = 2.5, base_height = 2.3)


####stat report####
df_t<-sumdf_long[which(sumdf_long$t_v=="training"),]
df_v<-sumdf_long[which(sumdf_long$t_v=="validation"),]

genstat<-function(x,y,stat,outdir){
  xvec<-x[which(x$attribute==stat),6]
  yvec<-y[which(y$attribute==stat),6]
  ks<-ks.test(xvec,yvec)
  shapirox<-shapiro.test(xvec)
  shapiroy<-shapiro.test(yvec)
  t_equal<-t.test(xvec,yvec,alternative = "two.sided", var.equal = TRUE)
  t_unequal<-t.test(xvec,yvec,alternative = "two.sided", var.equal = FALSE)
  wilcox<-wilcox.test(xvec,yvec,alternative = "two.sided")
  returnlist<-list(paste0("X is ", unique(x$t_v)),
                   paste0("Y is ", unique(y$t_v)),
                   paste0("Statistic is ",stat),
                   shapirox,shapiroy,ks,t_equal,t_unequal,wilcox)
  sink(paste0(outdir,"/stat_report_",stat,".txt"))
  return(returnlist)
}

genstat(x = df_t, y = df_v, stat = "variants_over_1%", outdir = "output_validation");sink()
genstat(x = df_t, y = df_v, stat = "variants_over_2%", outdir = "output_validation");sink()
genstat(x = df_t, y = df_v, stat = "Shannon_Entropy", outdir = "output_validation");sink()
genstat(x = df_t, y = df_v, stat = "Simpson_Diversity", outdir = "output_validation");sink()

#stop cluster
stopCluster(cl)

####clean up folders####
unlink("temp", recursive = TRUE)
unlink("fastq", recursive = TRUE)
dir.create("fastq")
