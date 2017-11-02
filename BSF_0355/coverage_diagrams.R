#setwd("/Volumes//vetgrid01/LCMV_project/Run_0277/Mappings/Coverages")
#setwd("~/Data/")
# test if there is exactly one argument and whether it is a file to be read
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1){
  write("An input filename containing coverage data must be supplied.\n",stderr())
  quit("no",status=1)
}
infile = args[1]

if (!file.exists(infile)) {
  write(paste("file",infile,"does not exist.\n"),stderr())
  quit("no",status=1)
}

virus.covs=read.table(infile,header=T)
outdir=paste0(dirname(infile),"/")

l.covs=virus.covs[virus.covs$CHROM=="L",2:ncol(virus.covs)]
s.covs=virus.covs[virus.covs$CHROM=="S",2:ncol(virus.covs)]
l.covs$windows=cut(l.covs$BPS,seq(0,max(l.covs$BPS),by=50))
s.covs$windows=cut(s.covs$BPS,seq(0,max(s.covs$BPS),by=50))

l.covs.50bp=aggregate(l.covs[,1:(ncol(l.covs)-1)],list(l.covs$win),mean)
s.covs.50bp=aggregate(s.covs[,1:(ncol(s.covs)-1)],list(s.covs$win),mean)


library(ggplot2)
library(RColorBrewer)
library("reshape2")
l.covs.50bp.melted <- melt(l.covs.50bp[,c(2,3:ncol(l.covs.50bp))],id.vars="BPS",variable.name = "sample", value.name="coverage")
s.covs.50bp.melted <- melt(s.covs.50bp[,c(2,3:ncol(s.covs.50bp))],id.vars="BPS",variable.name = "sample", value.name="coverage")

# too many samples, take n at a time
n=10
samples = sort(levels(l.covs.50bp.melted$sample))
s.plot = split(samples, floor((rank(samples)-1)/n))

for (i in seq(1,length(s.plot))) {
  bpL = ggplot(data=l.covs.50bp.melted[l.covs.50bp.melted$sample %in% s.plot[[i]],], aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
    geom_line() + scale_y_log10(breaks=c(10,50,100,250,500,1000,5000,10000,100000,500000)) + scale_linetype_manual(values = c(rep("dashed", 5), rep("solid", 5))) +
    scale_color_manual(values = c("black",brewer.pal(9, "Set1"))) + theme_bw() + 
    guides(linetype=guide_legend(ncol =2)) + ggtitle("L Segment") 
  bpL + ggsave(paste0(outdir,"lsegment_",i,".pdf"),width=8,height=6)
  bpS = ggplot(data=s.covs.50bp.melted[s.covs.50bp.melted$sample %in% s.plot[[i]],], aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
    geom_line() + scale_y_log10(breaks=c(10,50,100,250,500,1000,5000,10000,100000,500000)) + scale_linetype_manual(values = c(rep("solid", 5), rep("dashed", 5))) +
    scale_color_manual(values = c("black",brewer.pal(9, "Set1"))) + theme_bw() + 
    guides(linetype=guide_legend(ncol=2)) + ggtitle("S Segment") + ggsave("ssegment_low_covs.pdf",width=8,height=6)
  bpS + ggsave(paste0(outdir,"ssegment_",i,".pdf"),width=8,height=6)
}

q("no")

# garbage to be fitted in later
# get samples greater than 50K
l.covs.50bp.melted$sample[l.covs.50bp.melted$coverage > 50000]
unique(s.covs.50bp.melted$sample[s.covs.50bp.melted$coverage > 100000])
S31 S37 S38 S39
