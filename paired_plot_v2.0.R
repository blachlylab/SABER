library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr)
library(sqldf)
library(CrispRVariants) 
library(sangerseqR)
library("Rsamtools")
library(ggplot2)
library(rtracklayer)
library(GenomicFeatures)
library(gridExtra)
library(rtracklayer)
library(GenomicRanges)
library(RColorBrewer)
library(plyr)
library(rPython)
library(FField)
library(ggrepel)
library(cowplot)

#assign the samples
noumi_sample<-"AK005"
yesumi_sample<-"AK010"

# load the csv files


noumi<-read.table(paste0("output_noumi/thresh_0_tvaf_100_tpaf_0.1/no_UMI_",noumi_sample,".csv"), header=TRUE, sep=",")
yesumi<-read.table(paste0("output_yesumi/thresh_0_tvaf_100_tpaf_0.1/yes_UMI_",yesumi_sample,".csv"), header=TRUE, sep=",")

df1<-noumi
df2<-yesumi

#drop reads with edits at position -16 because they are junk

df1<-filter(df1, !grepl("-16", X))
df2<-filter(df2, !grepl("-16", X))

#keep the first N variants and drop the columns we don't need 
N<-22
df1.1<-df1[1:N,c(2,3,4,5,7,8)]
df2.1<-df2[1:N,c(2,3,4,5,7,8)]
colnames(df1.1)<-c("reads", "rnames", "is_common_variant", "is_no_variant", "specimen", "index")
colnames(df2.1)<-c("reads", "rnames", "is_common_variant", "is_no_variant", "specimen", "index")

#recalculate the allele frequencies without dropped alleles
df1.1$allele_freq<-df1.1$reads/sum(df1.1$reads)
df2.1$allele_freq<-df2.1$reads/sum(df2.1$reads)

# bind the dfs
df3<-full_join(df1.1, df2.1, by = "rnames")


#make palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}  
colors<-as.factor(gg_color_hue(13))
colors<-sort(colors)
colors

#apply colors to df
df3$color<-rep_len(colors, length.out = nrow(df3))
df3.1<-df3

#make unmatched variants white
levels<-levels(df3.1$color)
levels[length(levels)+1]<-"#000000"
df3.1$color<-factor(df3.1$color, levels=levels)
df3.1$color[is.na(df3.1$specimen.x)|is.na(df3.1$specimen.y)]<-"#000000"

#tidy up the data for plotting
df3.2<-df3.1[c(2,5,7,11,13,14)]
df3.3<-gather(data = df3.2, key = "attribute", value = "value_to_plot", allele_freq.x,allele_freq.y)

#plot the data

p3<-ggplot(df3.3, aes(x=attribute, y=value_to_plot, group = rnames)) +
  geom_point(aes(color=color, fill=color),size=2.5, alpha=0.4, shape=21, position=position_dodge(width=0.1)) +
  geom_line(aes(color=color),size=0.5, alpha=0.4, position=position_dodge(width=0.1)) +
  scale_y_log10()+
  scale_colour_identity()+
  scale_fill_identity()+
  scale_x_discrete(labels=c("Standard PCR","UMI PCR"))+
  theme_cowplot(font_size = 11)+
  labs(x="", y="Variant Allele Frequency")
p3
save_plot(plot = p3, filename = paste0("sidebyside_",noumi_sample,"_",yesumi_sample,"_7c6e3.pdf"), base_width = 2.5, base_height = 2.3)

#generate linear model
fit<-lm(allele_freq.y ~ allele_freq.x, data = df3.2)
summary(fit)
lb1<-round(summary(fit)$adj.r.squared,3)

#plot linerar regression
p1<-ggplot(df3.2, aes(x=allele_freq.x, y=allele_freq.y))
p1<-p1+
  geom_point(aes(color=color, fill=color), size=2.5, alpha=0.4,shape=21)+
  scale_fill_identity()+
  scale_color_identity()+
  geom_abline(intercept = fit[1]$coefficients[1], slope=fit[1]$coefficients[2], color="red", linetype="dashed")+
  annotate("text",x=0.3,y=0.05,label=parse(text = paste0('"Adjusted"','~R^2',":  ",lb1)))+
  theme_cowplot(font_size = 11)+
  labs(x="Standard PCR VAF", y="UMI PCR VAF")
p1
save_plot(plot = p1, filename = paste0("regression_",noumi_sample,"_",yesumi_sample,"_4597c.pdf"), base_width = 2.5, base_height = 2.3)

#generate data for qqplot
vec<-residuals(fit)
y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
x <- qnorm(c(0.25, 0.75))
slope <- diff(y)/diff(x)
int <- y[1L] - slope * x[1L]

d <- data.frame(resids = vec)
rnames_d<-rownames(d)
d<-cbind(rnames_d,d)

g<-ggplot(d, aes(sample = resids))+stat_qq()
d1<-ggplot_build(g)$data[[1]]
d1$rnames<-d$rnames_d[order(d$resids)]


levels(d1$rnames)
d1<-d1[order(d1$rnames),]

qq<-ggplot(d1,aes(theoretical,sample))
qq<-qq+
  geom_point(color = "black", fill = "white", size=2.5, alpha=0.4,shape=21)+
  geom_abline(slope = slope, intercept = int, color="black", linetype="dotted")+
  scale_color_identity()+
  scale_fill_identity()+
  theme_cowplot(font_size = 11)
qq
save_plot(plot = qq, paste0("qq_",noumi_sample,"_",yesumi_sample,"_78725.pdf"), base_width = 2.5, base_height = 2.3)



