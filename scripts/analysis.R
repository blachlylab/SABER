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
library(jsonlite)


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
# analysis_v2.0 cleans up analysis_v1.0

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
#enter the experiment data

genome<-"gestalt_pipeline4.fasta" # include .fasta.  Genome file has to be in genome folder.
genome_fp<-paste0(getwd(),"/references/",genome)
output_file<-"training_set"#for the cripsrvariants plot and prefix on sample output csv and vaf plot
group_name<-"training_set"

threshold<-200# number of reads below which they don't appear on the big crispr plot

#use_UMI<-FALSE#works how you think it should work
use_UMI<-fromJSON("config.json")$use_umi


####core pipeline####
#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


# make working subdirectories and create variables 
output_folder<-"output"
fastqc_files_folder<-paste0(output_folder,"/fastqc_files")
bam_output_folder<-paste0(output_folder,"/bam")
PEAR_output_folder<-paste0(output_folder,"/pear")
fastq_trimmed<-paste0(getwd(),"/fastq_trimmed")
fastq_merged<-paste0(getwd(),"/fastq_merged")
cutadapt1_folder<-paste0(getwd(),"/cutadapt5p")
cutadapt2_folder<-paste0(getwd(),"/cutadapt3p")
fastq_umi_out<-paste0(getwd(),"/fastq_umi")
bam_dedup<-paste0(getwd(),"/dedup")
bam_dedup_output_folder<-paste0(output_folder,"/dedup")

bam_files<-list.files(path = bam_output_folder,".*.bam$"); bam_fp<-paste0(bam_output_folder,"/",bam_files)

bam_files

mid_names<-substr(bam_files, 6,10)

if (use_UMI==TRUE){
    bam_dedup_files<-list.files(path = bam_dedup_output_folder,".*.bam$"); bam_dedup_fp<-paste0(bam_dedup_output_folder,"/",bam_dedup_files)
}

#build metadata for experiment
if (use_UMI==TRUE){
  bam_fnames<-bam_dedup_fp
} else {
  bam_fnames <- bam_fp
}
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

#while (!is.null(dev.list())) dev.off()

# TODO:configurable height/width?
#dev.copy2pdf(file=paste0(output_folder,"/",output_file,".pdf"), width = 36, height = 36)  #for 10 samples use 24 x 24, for 30-40 samples use 48 x 48
pdf(file=paste0(output_folder,"/",output_file,".pdf"), width = 36, height = 36)
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
dev.off()






#####generate informative, no variant and common variant tables with vaf and paf thresholds####
generate_v1<-function(x,thresh,vaf_thresh,paf_thresh,vc_prop,output_folder,output_file){
  df<-x#isolate df 
  rnames<-row.names(df)#duplicate row names as new variable
  df$rnames<-rnames#add row name variabe to df
  df<-df[which(df$rnames!="Other"),]# removes other
  ##important:  CrispRVariants proportion generates PERCENT values, not FRACTION values.  So TVAF = 1 => 1% vaf, not 100% VAF.
  common_vars<-names(which(rowSums(vc_prop>vaf_thresh)/ncol(vc_prop)>(paf_thresh)))#returns the names of rows with a propotion greater than vaf_thresh occuring in more than paf_thresh
  ##end important
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
vc_prop<-variantCounts(crispr_set, result = "proportions")# generates a matrix of variant allele percents, similar to vc_all


thresh<-0;tvaf<-100;tpaf<-0.05;parLapply(cl,vc_list,generate_v1,thresh = thresh, vaf_thresh = tvaf, paf_thresh = tpaf,vc_prop = vc_prop,output_folder = output_folder, output_file = output_file)
thresh<-0;tvaf<-10;tpaf<-0.05;parLapply(cl,vc_list,generate_v1,thresh = thresh, vaf_thresh = tvaf, paf_thresh = tpaf,vc_prop = vc_prop,output_folder = output_folder, output_file = output_file)
thresh<-0;tvaf<-2;tpaf<-0.05;parLapply(cl,vc_list,generate_v1,thresh = thresh, vaf_thresh = tvaf, paf_thresh = tpaf,vc_prop = vc_prop,output_folder = output_folder, output_file = output_file)
thresh<-0;tvaf<-1;tpaf<-0.05;parLapply(cl,vc_list,generate_v1,thresh = thresh, vaf_thresh = tvaf, paf_thresh = tpaf,vc_prop = vc_prop,output_folder = output_folder, output_file = output_file)
thresh<-0;tvaf<-0.1;tpaf<-0.05;parLapply(cl,vc_list,generate_v1,thresh = thresh, vaf_thresh = tvaf, paf_thresh = tpaf,vc_prop = vc_prop,output_folder = output_folder, output_file = output_file)
thresh<-0;tvaf<-100;tpaf<-0.1;parLapply(cl,vc_list,generate_v1,thresh = thresh, vaf_thresh = tvaf, paf_thresh = tpaf,vc_prop = vc_prop,output_folder = output_folder, output_file = output_file)
thresh<-0;tvaf<-10;tpaf<-0.1;parLapply(cl,vc_list,generate_v1,thresh = thresh, vaf_thresh = tvaf, paf_thresh = tpaf,vc_prop = vc_prop,output_folder = output_folder, output_file = output_file)
thresh<-0;tvaf<-2;tpaf<-0.1;parLapply(cl,vc_list,generate_v1,thresh = thresh, vaf_thresh = tvaf, paf_thresh = tpaf,vc_prop = vc_prop,output_folder = output_folder, output_file = output_file)
thresh<-0;tvaf<-1;tpaf<-0.1;parLapply(cl,vc_list,generate_v1,thresh = thresh, vaf_thresh = tvaf, paf_thresh = tpaf,vc_prop = vc_prop,output_folder = output_folder, output_file = output_file)
thresh<-0;tvaf<-0.1;tpaf<-0.1;parLapply(cl,vc_list,generate_v1,thresh = thresh, vaf_thresh = tvaf, paf_thresh = tpaf,vc_prop = vc_prop,output_folder = output_folder, output_file = output_file)

#stop cluster
stopCluster(cl)


####sparsity####
#calculate sparsity of the original dataframe
#load the data
sparsity<-as.data.frame(colSums(vc_all==0)/(colSums(vc_all==0) + colSums(vc_all!=0)))#95-99% sparse
colnames(sparsity)<-"Sparsity"
sparsity<-sparsity %>%
  rownames_to_column("Sample")
sparsity$nonsparse<-1-sparsity$Sparsity
sparsity$index<-factor(rank(sparsity$nonsparse, ties.method = "first"))

# plot sparsity
sparse_plot<-ggplot(sparsity, aes(x = sparsity$index, y = sparsity$Sparsity))
sparse_plot<-sparse_plot+
  geom_bar(stat = "identity", color = "black", fill = "#fde725")+
  scale_x_discrete(breaks = c(5,10,15,20,25,30,35))+
  theme_cowplot(font_size = 11)+
  ylab(label = "Sparsity")+
  xlab(label = "Sample Index")
  #theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
sparse_plot
save_plot(plot = sparse_plot, filename = "output/sparsity_plot_e7ba6.pdf", base_height = 2, base_width = 3.5)


####big heatmap####
# plot heatmap of all variants present over 20000 reads in any sample.  i.e. modify the crispr plots heatmap
heat1_df<-vc_all
heat1_df<-heat1_df %>% 
  rownames_to_column("Variants") %>%
  filter_if(is.numeric, any_vars(. >20000))%>%
  column_to_rownames("Variants")

#cluster
heat1_matrix<-as.matrix(heat1_df)
heat1_dendro<-as.dendrogram(hclust(d=dist(x=heat1_matrix)))

x <- heat1_matrix
dd.col <- as.dendrogram(hclust(dist(x)))
dd.row <- as.dendrogram(hclust(dist(t(x))))
dx <- dendro_data(dd.row)
dy <- dendro_data(dd.col)

# helper function for creating dendograms
ggdend <- function(df) {
  ggplot() +
    geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
    labs(x = "", y = "") + theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank())
}

# x/y dendograms
px <- ggdend(dx$segments)
py <- ggdend(dy$segments) + coord_flip()

#px
#py

heat1_df<-heat1_df %>%
  rownames_to_column("Variants")
heat1_df_long<-melt(heat1_df)

# #reorder heatmap
heat1_orderx<-order.dendrogram(dd.col)
heat1_ordery<-order.dendrogram(dd.row)
heat1_df_long$Variants<-factor(x=heat1_df_long$Variants, levels = heat1_df$Variants[rev(heat1_orderx)], ordered = TRUE)
heat1_df_long$variable<-factor(x=heat1_df_long$variable, levels = colnames(vc_all)[heat1_ordery], ordered = TRUE) 
heat1_df_long$rank<-rank(heat1_df_long$variable, ties.method = "min")
heat2_df_long<-unique(heat1_df_long[,c(2,4)])
heat2_df_long<-arrange(heat2_df_long,rank)
heat2_df_long$index<-factor(seq.int(nrow(heat2_df_long)))
heat3_df_long<-left_join(x = heat1_df_long, y = heat2_df_long, by = "variable")


heat1<-ggplot(data = heat3_df_long, aes(x = index, y=Variants, fill = value))
heat1<-heat1+
  geom_tile()+
  scale_fill_viridis_c(begin = 0, end = 1, guide = "colourbar", aesthetics = "fill", limits = c(0,60000), option = "D")+
  scale_x_discrete(breaks = c(5,10,15,20,25,30,35))+
  theme_cowplot(font_size = 11)+ 
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y = element_blank())+
  theme(axis.ticks.y = element_blank())+
  theme(axis.line.y = element_blank())+
  guides(fill=guide_legend(title="Read Counts", reverse = T))+
  xlab("Sample Index")+
  ylab("Variant")+
  coord_fixed()
heat1
save_plot(plot = heat1, filename = "output/full_heatmap_b830e.pdf", base_width = 4, base_aspect_ratio = 1.3)##check to be sure what the cutoff is!!

#####analysis of tvaf, cvsa settings####
cor_func2<-function(include_uninformative, upper_lower, thresh, tvaf, tpaf, method, tag){
  input_directory<-paste0("output/thresh_",thresh,"_tvaf_",tvaf,"_tpaf_",tpaf)
  inlist<-list.files(input_directory,full.names = TRUE, pattern = "\\d.csv$")
  df_cor<-data.frame(Variant = character(), Count = integer(), Specimen = character())
  for (i in 1:length(inlist)){
    df_new<-read.csv(file = inlist[i], header = TRUE)
    df_new<-df_new[,c(3,2,7)]
    colnames(df_new)<-c("Variant", "Count", "Specimen")
    df_cor<-rbind(df_cor,df_new)
  }
  df_cor_wide<-reshape(df_cor, idvar = "Variant", timevar = "Specimen", direction = "wide")
  df_cor_wide[is.na(df_cor_wide)]<-0
  colnames(df_cor_wide)<-substr(colnames(df_cor_wide),7,11)
  rownames(df_cor_wide)<-df_cor_wide[,1]
  df_cor_wide<-df_cor_wide[,-1]
  if (include_uninformative==TRUE) {df_cor_wide1<-df_cor_wide} else {df_cor_wide1<-df_cor_wide[-c(1,2),]}
  cormat1<-round(cor(df_cor_wide1, method = method),2)
  #p<-cor.test(df_cor_wide1, method = method)
  #helper functions
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(x){
    x[upper.tri(x)] <- NA
    return(x)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(x){
    x[lower.tri(x)]<- NA
    return(x)
  }
  reorder_cormat <- function(x){
    # Use correlation between variables as distance
    dd <- as.dist((1-x)/2)
    hc <- hclust(dd)
    x <-x[hc$order, hc$order]
  }
  # Reorder the correlation matrix
  cormat2 <- reorder_cormat(cormat1)
  cormat3 <- as.data.frame(cormat2)
  colnames(cormat3)<-c(1:ncol(cormat3));rownames(cormat3)<-c(1:nrow(cormat3))
  upper_tri <- get_upper_tri(cormat2)
  lower_tri <-get_lower_tri(cormat2)
  upper_tri1<-get_upper_tri(cormat3)
  if (upper_lower=="upper") {melted_cormat<-melt(upper_tri, na.rm = TRUE)} else {melted_cormat<-melt(lower_tri, na.rm = TRUE)}
  if (upper_lower=="upper") {melted_cormat1<-melt(upper_tri1, na.rm = TRUE)} else {melted_cormat1<-melt(lower_tri1, na.rm = TRUE)}
  melted_cormat2<-cbind(melted_cormat,melted_cormat1)
  colnames(melted_cormat2)<-c("sample_y","sample_x","value1","index","value2")
  
  #calculate Fraction informative
  FI<-as.data.frame(colSums(df_cor_wide[-c(1,2),])/colSums(df_cor_wide))#ratio of reads excluding no variant and common variant to all reads in each sample.
  colnames(FI)<-"Fraction_Informative"
  FI$sample_x<-rownames(FI)
  
  #join dfs
  cordf<-left_join(x = melted_cormat2, y = FI, by = "sample_x")
  
    #plot the fi data - not using the plot currently
  fidf<-unique(cordf[,c(2,6)])
  fiplot<-ggplot(fidf, aes(x=fidf$sample_x, y = fidf$Fraction_Informative))
  fiplot<-fiplot+
    geom_bar(stat="identity", color = "black", fill = "black")+
    theme_cowplot(font_size = 8)+
    theme(axis.text.x = element_blank())+
    theme(axis.line.x = element_blank())+
    theme(axis.ticks.x = element_blank())+
    #theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    ylab(label = "Fraction\nInformative")+
    ylim(0,1)+
    xlab(label = "")
  fiplot
  
  #plot the data
  avg_cor<-round(mean(cordf[which(cordf$sample_y!=cordf$sample_x),3]),4)#averages pearson correlation excluding identical comparisons
  mean_inf<-round(mean(fidf$Fraction_Informative),3)
  if (method=="pearson") {corlabel<-"Pearson"}; if (method=="spearman") {corlabel<-"Spearman"}
  corplot2<-ggplot(data = cordf, aes(x = cordf$index, y = cordf$sample_y, fill = cordf$value1))
  corplot2<-corplot2+
    geom_tile()+
    scale_fill_viridis_c(begin = 0, end = 1, guide = "colourbar", aesthetics = "fill")+
    scale_x_discrete(breaks = c(5,10,15,20,25,30,35))+
    theme_cowplot(font_size = 11)+ 
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
    theme(axis.text.y = element_blank())+
    #theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    theme(axis.line.y = element_blank())+
    theme(axis.ticks.y = element_blank())+
    xlab("Sample Index")+
    ylab("")+
    guides(fill=guide_legend(title=paste0(corlabel,"\nCorrelation"), reverse = T))+
    theme(legend.position = c(0,0.7))+
    coord_fixed()
  corplot2
  corplot3<-plot_grid(fiplot,corplot2, align = "v", nrow = 2, rel_heights = c(1/4, 3/4))
  save_plot(plot = corplot2, filename = paste0("output/corplot2_thresh_",thresh,"_tvaf_",tvaf,"_tpaf_",tpaf,"_uninform_",include_uninformative,"_",tag,".pdf"), base_width = 2.5, base_height = 2.75)
  nubbin<-data.frame(threshold = thresh, tvaf = tvaf, tpaf = tpaf, mean_inf = mean_inf, avg_cor = avg_cor)
  write.csv(x = nubbin, file = paste0(input_directory,"/nubbin_",include_uninformative,"_",method,".csv"))
}

cor_func2(include_uninformative = TRUE,upper_lower = "upper", thresh = 0, tvaf = 100, tpaf = 0.05, method = "pearson", tag = "2cbda")


cor_func2(include_uninformative = FALSE,upper_lower = "upper", thresh = 0, tvaf = 0.1, tpaf = 0.05, method = "pearson", tag = "4gD9x")
cor_func2(include_uninformative = FALSE,upper_lower = "upper", thresh = 0, tvaf = 1, tpaf = 0.05, method = "pearson", tag = "837p0")
cor_func2(include_uninformative = FALSE,upper_lower = "upper", thresh = 0, tvaf = 2, tpaf = 0.05, method = "pearson", tag = "89eo6")
cor_func2(include_uninformative = FALSE,upper_lower = "upper", thresh = 0, tvaf = 10, tpaf = 0.05, method = "pearson", tag = "3qzID")
cor_func2(include_uninformative = FALSE,upper_lower = "upper", thresh = 0, tvaf = 100, tpaf = 0.05, method = "pearson", tag = "124fD")
cor_func2(include_uninformative = FALSE,upper_lower = "upper", thresh = 0, tvaf = 0.1, tpaf = 0.1, method = "pearson", tag = "9Egr0")
cor_func2(include_uninformative = FALSE,upper_lower = "upper", thresh = 0, tvaf = 1, tpaf = 0.1, method = "pearson", tag = "GMvgC")
cor_func2(include_uninformative = FALSE,upper_lower = "upper", thresh = 0, tvaf = 2, tpaf = 0.1, method = "pearson", tag = "0G6O6")
cor_func2(include_uninformative = FALSE,upper_lower = "upper", thresh = 0, tvaf = 10, tpaf = 0.1, method = "pearson", tag = "Q7Fv4")

####plot correlation vs fraction informative####
pearson_inlist<-list.files("output", recursive = TRUE, pattern = "nubbin_FALSE_pearson.csv",full.names = TRUE)
pearson_inlist1<-as.character(sort(factor(pearson_inlist, levels = factor(pearson_inlist[c(1:4,8,9,5:7)]))))

curve_df0<-read.csv(pearson_inlist1[1])
for (i in 2:length(pearson_inlist1)){
  curve_df0<-rbind(curve_df0,read.csv(pearson_inlist1[i]))
}

curve_df1<-curve_df0

curve_df2<-curve_df1 %>%
  rownames_to_column("condition")
curve_df2<-curve_df2[,-2]
curve_df2[,1] <- as.factor(curve_df2[,1])

legend_text<-function(x){paste0("VAF>",curve_df2[x,3]/100," & PAF>",curve_df2[x,4])}# this adjusts the VAF label to correspond to FRACTIONS instead of PERCENT.  Data unchanged.

curve_plot<-ggplot(curve_df2, aes(x = curve_df2$avg_cor, y = curve_df2$mean_inf))
curve_plot<-curve_plot+
  geom_jitter(shape = 21, size = 3, color = "black", alpha = 0.8, aes(fill = curve_df2$condition))+
  theme_cowplot(font_size = 11)+
  scale_fill_brewer(palette = "Set1", labels = c(legend_text(1),
                                                 legend_text(2),
                                                 legend_text(3),
                                                 legend_text(4),
                                                 legend_text(5),
                                                 legend_text(6),
                                                 legend_text(7),
                                                 legend_text(8),
                                                 legend_text(9)), name = "")+
  theme(legend.position = c(0.40,0.43))+
  labs(x = "Mean Correlation Between Samples", y = "Mean Fraction Informative", fill = "Condition")
curve_plot
save_plot(plot = curve_plot, filename = paste0(output_folder,"/curve_plot_8Ty02.pdf"), base_height = 3, base_aspect_ratio = 1.1)


####re-run training samples using condition tvaf = 1, tpaf = 0.1 and identify informative vs uninformative samples using kmeans####
#load the data
infiles<-list.files("output/thresh_0_tvaf_1_tpaf_0.1", pattern = "\\d.csv", full.names = TRUE)

fidf_0<-read.csv(infiles[1])
fi<-1-sum(c(fidf_0[1,2],fidf_0[2,2]))/sum(fidf_0[2])
fidf_1<-fidf_0[,c(7,2,3)]
colnames(fidf_1)<-c("specimen", "counts", "variant")
fidf_1$FI<-rep(x = fi,times = nrow(fidf_0))
for (i in 2:length(infiles)){
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


md<-read.csv(list.files("references", pattern = "clone_md.csv", full.names = TRUE))#reads in the metadata to identify sample groups
colnames(md)<-c("specimen", "class", "batch")
fidf_4<-merge(fidf_3,md,by = "specimen", all.x = TRUE)

fidf_3.1<-fidf_3

#plot the fi distribution with kmeans
fi_dist<-ggplot(fidf_3, aes(x=factor(fidf_3$rank), y=fidf_3$FI, fill=fidf_3$cluster))
fi_dist<-fi_dist+
  geom_bar(stat = "identity", color = "black")+
  scale_fill_manual(values = c("#72ce55","#fde725"),name = "")+
  scale_x_discrete(breaks = c(5,10,15,20,25,30,35))+
  theme_cowplot(font_size = 11)+
  theme(legend.position = c(0.6,0.95))+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  labs(x="Sample Index", y="Fraction Informative")
fi_dist
save_plot(plot = fi_dist, filename = "output/fidist_6e65e.pdf", base_width = 4, base_height = 3)

#batch plot
md<-read.csv(list.files("references", pattern = "clone_md.csv", full.names = TRUE))
colnames(md)<-c("specimen", "class", "batch","t_v")
fidf_4<-merge(fidf_3,md,by = "specimen", all.x = TRUE)
fidf_4$batch<-as.factor(fidf_4$batch)
fi_batch<-ggplot(fidf_4, aes(x=reorder(fidf_4$specimen,fidf_4$rank), y=fidf_4$FI, fill=fidf_4$batch))
fi_batch<-fi_batch+
  geom_bar(stat = "identity", color = "black")+
  scale_fill_manual(values = c("#3C5488","#DC0000"),
                    name = "",
                    breaks=c("1","2"),
                    labels=c("Batch 1","Batch 2"))+
  theme_cowplot(font_size = 11)+
  theme(legend.position = c(0.6,0.95))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8))+
  labs(x="Sample", y="Fraction Informative")
fi_batch
save_plot(plot = fi_batch, filename = paste0(output_folder,"/fibatch_6b3e2.pdf"), base_height = 3, base_aspect_ratio = 1.3)

all_samples<-fidf_4$specimen
informative_samples<-filter(fidf_4, fidf_4$cluster=="Informative")

####quantifying alleles from thresholded data in informative samples####
#reassemble the long form data frame of read counts using only informative samples
inlist<-as.data.frame(cbind(as.character(all_samples), infiles))
colnames(inlist)<-c("specimen", "path")
inlist_informative<-left_join(informative_samples,inlist, by = "specimen")
informative_path<-as.character(inlist_informative$path)
qdf0<-read.csv(informative_path[1])
colnames(qdf0)<-c("X","read_count","rnames","is_common_variant","is_no_variant","allele_freq","specimen","index")
for (i in 2:length(informative_path)) {
  qdf_new<-read.csv(informative_path[i])
  colnames(qdf_new)<-c("X","read_count","rnames","is_common_variant","is_no_variant","allele_freq","specimen","index")
  qdf0<-rbind(qdf0,qdf_new)
}
qdf1<-qdf0[,-1]

#generate ordered distribution of variant alleles and plot at various cutoffs
##all vars so no cutoff - this generates plots for all samples
dir.create("output/thresh_0_tvaf_1_tpaf_0.1/all_vars"); outdir1<-"output/thresh_0_tvaf_1_tpaf_0.1/all_vars"
spec_list<-unique(qdf1$specimen)
for (i in 1:length(spec_list)){
  qdf_filt<-filter(qdf1, qdf1[,6]==spec_list[i])#loads in dfs one sample at a time
  qdf_filt_inf<-qdf_filt[-c(1,2),]# removes no variant and common variant
  qdf_filt_inf$index<-qdf_filt_inf$index-2#resets index count
  p<-ggplot(data = qdf_filt_inf, aes(x = qdf_filt_inf$index, y = qdf_filt_inf$allele_freq))
  p<-p+
    theme_cowplot(font_size = 11)+
    geom_bar(stat = "identity", color = "black")+
    labs(x="index", y = "Raw VAF", title = paste0(qdf_filt_inf[1,6],"\nRichness=",nrow(qdf_filt_inf)))#does not recalculate VAF with exclusion of no variant and common variant
  save_plot(filename = paste0(outdir1,"/",qdf_filt_inf[1,6],"_all_vars.pdf"), plot = p, base_height = 3, base_aspect_ratio = 1.3)
}

# make a plot of read counts for a single example for paper
qdf_AB029<-filter(qdf1,qdf1$specimen=="AB029")
qdf_AB029$log10_read_count<-log10(qdf_AB029$read_count)
qdf_AB029_inf<-qdf_AB029[-c(1,2),];qdf_AB029_inf$index<-qdf_AB029_inf$index-3
qdf_AB029_inf_filt<-filter(qdf_AB029_inf, qdf_AB029_inf$index<25)
AB029_plot<-ggplot(data = qdf_AB029_inf_filt, aes(x = qdf_AB029_inf_filt$index, y = qdf_AB029_inf_filt$log10_read_count))
AB029_plot<-AB029_plot+
  theme_cowplot(font_size = 11)+
  geom_bar(stat = "identity", color = "black", fill = "#DC0000")+
  labs(x="Variant Index", y = expression(paste(Log[10]," Read Count",sep="")))#does not recalculate VAF with exclusion of no variant and common variant
AB029_plot
save_plot(filename = "output/AB029_example_758ed.pdf", plot = AB029_plot, base_width = 2.5, base_height = 2.3)

#generate a vector of shannon entropy values for the informative training set
shannon<-numeric()
for (i in 1:length(spec_list)){
  qdf_filt<-filter(qdf1, qdf1[,6]==spec_list[i])#loads in dfs one sample at a time
  qdf_filt_inf<-qdf_filt[-c(1,2),]# removes no variant and common variant
  qdf_filt_inf$index<-qdf_filt_inf$index-2#resets index count
  shannon<-c(shannon,shannon_func(qdf_filt_inf[1]))#calculates off of read count
}
shannon
shapiro.test(shannon)

#generate a vector of simpson diversity indices for the informative training set
simpson<-numeric()
for (i in 1:length(spec_list)){
  qdf_filt<-filter(qdf1, qdf1[,6]==spec_list[i])#loads in dfs one sample at a time
  qdf_filt_inf<-qdf_filt[-c(1,2),]# removes no variant and common variant
  qdf_filt_inf$index<-qdf_filt_inf$index-2#resets index count
  simpson<-c(simpson,simpson_func(qdf_filt_inf[1]))#calculates off of read count
}
simpson
shapiro.test(simpson)

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
shapiro.test(vaf2)

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
shapiro.test(vaf1)

#compile and summarize data
sumdf<-as.data.frame(cbind(levels(qdf1$specimen),vaf1,vaf2,shannon,simpson))
colnames(sumdf)<-c("Sample", "variants_over_1%", "variants_over_2%", "Shannon_Entropy", "Simpson_Diversity")

sumdf_long<-gather(sumdf, key = "attribute", value = "value_to_plot", -Sample)
sumdf_long$set<-c(rep("training_set", times = nrow(sumdf_long)))
sumdf_long$value_to_plot<-as.numeric(sumdf_long$value_to_plot)
vaf1_long<-filter(sumdf_long, sumdf_long$attribute=="variants_over_1%")
vaf2_long<-filter(sumdf_long, sumdf_long$attribute=="variants_over_2%")
shannon_long<-filter(sumdf_long, sumdf_long$attribute=="Shannon_Entropy")
simpson_long<-filter(sumdf_long, sumdf_long$attribute=="Simpson_Diversity")


p1<-ggplot(vaf1_long, aes(x = vaf1_long$value_to_plot))
p1<-p1+
  geom_histogram(aes(y = ..density..), color = "black", fill = "#DC0000", binwidth = 1) + 
  geom_density(alpha=.2, fill="#DC0000")+
  theme_cowplot(font_size = 11)+
  theme(legend.position = "none")+
  scale_x_continuous(limits = c(0,18), breaks = c(0,2,4,6,8,10,12,14,16,18))+
  labs(x="Clones with VAF>0.01", y="Density")
p1
save_plot(plot = p1, filename = "output/vaf_0.01_646d9.pdf", base_width = 2.5, base_height = 2.3)


p2<-ggplot(vaf2_long, aes(x = vaf2_long$value_to_plot))
p2<-p2+
  geom_histogram(aes(y = ..density..), color = "black", fill = "#DC0000", binwidth = 1) + 
  geom_density(alpha=.2, fill="#DC0000")+
  theme_cowplot(font_size = 11)+
  theme(legend.position = "none")+
  scale_x_continuous(limits = c(0,10.5), breaks = c(0,2,4,6,8,10))+
  labs(x="Clones with VAF>0.02", y="Density")
p2
save_plot(plot = p2, filename = "output/vaf_0.02_f0c32.pdf", base_width = 2.5, base_height = 2.3)

p3<-ggplot(shannon_long, aes(x = shannon_long$value_to_plot))
p3<-p3+
  geom_histogram(aes(y = ..density..), color = "black", fill = "#DC0000", binwidth = 0.5) + 
  geom_density(alpha=.2, fill="#DC0000")+
  theme_cowplot(font_size = 11)+
  theme(legend.position = "none")+
  #scale_x_continuous(limits = c(0,5.5), breaks = c(0,1,2,3,4,5))+
  labs(x="Shannon Entropy", y="Density")
p3
save_plot(plot = p3, filename = "output/shannon_95119.pdf", base_width = 2.5, base_height = 2.3)

p4<-ggplot(simpson_long, aes(x = simpson_long$value_to_plot))
p4<-p4+
  geom_histogram(aes(y = ..density..), color = "black", fill = "#DC0000", binwidth = 1) + 
  geom_density(alpha=.2, fill="#DC0000")+
  theme_cowplot(font_size = 11)+
  theme(legend.position = "none")+
  #scale_x_continuous(limits = c(0,7.5), breaks = c(0,1,2,3,4,5,6,7))+
  labs(x="Simpson Diversity", y="Density")
p4
save_plot(plot = p4, filename = "output/simpson_9780d.pdf", base_width = 2.5, base_height = 2.3)

####go back to whole training set (informative + uninformative) and plot fraction informative vs clones over 2%
inlist_informative2<-left_join(fidf_4[,-7],inlist, by = "specimen")#all samples included
informative_path2<-as.character(inlist_informative2$path)
qdf2<-read.csv(informative_path2[1])
colnames(qdf2)<-c("X","read_count","rnames","is_common_variant","is_no_variant","allele_freq","specimen","index")
for (i in 2:length(informative_path2)) {
  qdf_new<-read.csv(informative_path2[i])
  colnames(qdf_new)<-c("X","read_count","rnames","is_common_variant","is_no_variant","allele_freq","specimen","index")
  qdf2<-rbind(qdf2,qdf_new)
}
qdf3<-qdf2[,-1]

#generate a vector of the number of variants over 2%
vaf2_all<-numeric()
for (i in 1:length(all_samples)){
  qdf_filt<-filter(qdf3, qdf3[,6]==all_samples[i])#loads in dfs one sample at a time
  qdf_filt_inf<-qdf_filt[-c(1,2),]# removes no variant and common variant
  qdf_filt_inf$index<-qdf_filt_inf$index-2#resets index count
  qdf_filt_inf$adj_allele_freq<-qdf_filt_inf$read_count/sum(qdf_filt_inf$read_count)#recalculates allele frequencies from read counts without no variant and common variant
  qdf_2percent<-filter(qdf_filt_inf, qdf_filt_inf$adj_allele_freq>=0.02)
  vaf2_all<-c(vaf2_all,nrow(qdf_2percent))
}
vaf2_all
fidf_5<-fidf_4
fidf_5$vars_over_0.02<-vaf2_all

colorder<-c("Informative","Not Informative")
p5<-ggplot(fidf_5, aes(x = fidf_5$cluster, y = fidf_5$vars_over_0.02, fill = fidf_5$cluster))
p5<-p5+
  geom_boxplot()+
  theme_cowplot(font_size = 11)+
  scale_x_discrete(limits = colorder, labels = c("Informative", "Not\nInformative"))+
  scale_fill_manual(values = c("#72ce55","#fde725"))+
  theme(legend.position = "none")+
  labs(x = "", y = "Clones with VAF>0.02")
p5
save_plot(plot = p5, filename = "output/informative_boxplot_a1332.pdf", base_width = 2.5, base_height = 2.3)

# clean up directories
#dir.create("fastq")


