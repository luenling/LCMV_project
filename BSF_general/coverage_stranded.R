#setwd("/Volumes//vetgrid01/LCMV_project/Run_0277/Mappings/Coverages")
#setwd("~/Data/")
# test if there is exactly one argument and whether it is a file to be read
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2){
  write("2 input filenames containing coverage data for the plus and minus strand must be supplied.\n",stderr())
  quit("no",status=1)
}
pfile = args[1]
mfile = args[2]

if (!file.exists(pfile) | !file.exists(mfile)) {
  write(paste("file",pfile,"or", mfile, "does not exist.\n"),stderr())
  quit("no",status=1)
}

#pfile="/Volumes/vetlinux01/LCMV/Run_0355/BQSR/Coverages/bsf_0355_all_plus.depths"
#mfile="/Volumes/vetlinux01/LCMV/Run_0355/BQSR/Coverages/bsf_0355_all_minus.depths"

virus.p=read.table(pfile,header=T)
virus.m=read.table(mfile,header=T)
outdir=paste0(dirname(pfile),"/")

l.p=virus.p[virus.p$CHROM=="L",2:ncol(virus.p)]
s.p=virus.p[virus.p$CHROM=="S",2:ncol(virus.p)]
l.p$windows=cut(l.p$BPS,seq(0,max(l.p$BPS),by=50))
s.p$windows=cut(s.p$BPS,seq(0,max(s.p$BPS),by=50))
l.m=virus.m[virus.m$CHROM=="L",2:ncol(virus.m)]
s.m=virus.m[virus.m$CHROM=="S",2:ncol(virus.m)]
l.m$windows=cut(l.m$BPS,seq(0,max(l.m$BPS),by=50))
s.m$windows=cut(s.m$BPS,seq(0,max(s.m$BPS),by=50))

l.p.50bp=aggregate(l.p[,1:(ncol(l.p)-1)],list(l.p$win),mean)
s.p.50bp=aggregate(s.p[,1:(ncol(s.p)-1)],list(s.p$win),mean)
l.m.50bp=aggregate(l.m[,1:(ncol(l.m)-1)],list(l.m$win),mean)
s.m.50bp=aggregate(s.m[,1:(ncol(s.m)-1)],list(s.m$win),mean)

library(ggplot2)
library(RColorBrewer)
library("reshape2")
l.p.50bp.melted <- melt(l.p.50bp[,c(2,3:ncol(l.p.50bp))],id.vars="BPS",variable.name = "sample", value.name="coverage")
s.p.50bp.melted <- melt(s.p.50bp[,c(2,3:ncol(s.p.50bp))],id.vars="BPS",variable.name = "sample", value.name="coverage")
l.m.50bp.melted <- melt(l.m.50bp[,c(2,3:ncol(l.m.50bp))],id.vars="BPS",variable.name = "sample", value.name="coverage")
s.m.50bp.melted <- melt(s.m.50bp[,c(2,3:ncol(s.m.50bp))],id.vars="BPS",variable.name = "sample", value.name="coverage")
l.p.50bp.melted$strand="sense"
s.p.50bp.melted$strand="sense"
l.m.50bp.melted$strand="anti-sense"
s.m.50bp.melted$strand="anti-sense"
l.50bp.melted=rbind(l.p.50bp.melted,l.m.50bp.melted)
s.50bp.melted=rbind(s.p.50bp.melted,s.m.50bp.melted)
l.50bp.melted$strand=factor(l.50bp.melted$strand,levels=c("sense","anti-sense"))
s.50bp.melted$strand=factor(s.50bp.melted$strand,levels=c("sense","anti-sense"))

# too many samples, take n at a time
n=8
samples = sort(levels(l.p.50bp.melted$sample))
s.plot = split(samples, floor((rank(samples)-1)/n))

for (i in seq(1,length(s.plot))) {
  bpL = ggplot(data=l.50bp.melted[l.50bp.melted$sample %in% s.plot[[i]],], aes(x=BPS, y=coverage, group=strand, colour = strand)) +
    geom_line() + #scale_y_log10(breaks=c(10,50,100,250,500,1000,5000,10000,100000,500000))
    scale_color_manual(values=c("blue","red")) + theme_bw() +  facet_grid(sample ~ ., scales = "free") +
   ggtitle("L Segment") + theme(legend.position = "bottom",legend.title=element_blank())
  bpL + ggsave(paste0(outdir,"lsegment_stranded_",paste(s.plot[[i]],sep="_",collapse="_"),".pdf"),width=8,height=8)
  bpS = ggplot(data=s.50bp.melted[s.50bp.melted$sample %in% s.plot[[i]],], aes(x=BPS, y=coverage, group=strand, colour = strand)) +
    geom_line() + #scale_y_log10(breaks=c(10,50,100,250,500,1000,5000,10000,100000,500000))
    scale_color_manual(values=c("blue","red")) + theme_bw() +  facet_grid(sample ~ ., scales = "free") +
    ggtitle("S Segment") + theme(legend.position = "bottom",legend.title=element_blank())
  bpS + ggsave(paste0(outdir,"ssegment_stranded_",paste(s.plot[[i]],sep="_",collapse="_"),".pdf"),width=8,height=8)
}

q("no")
