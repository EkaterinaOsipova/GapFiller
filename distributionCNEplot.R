#!/usr/bin/env Rscript

######
setwd('/Users/osipova/Documents/R_scripts/')

gainedCNEs = read.table('gained.CNEs.G.bed', sep=' ', header=F)
lostCNEs = read.table('lost.CNEs.L.bed', sep=' ', header=F)
#transposons = read.table('allTransposons.hg38.bed', sep='\t', header=F)
simpleRepeats = read.table('LTRtransposons.bed', sep='\t', header=F)
total = rbind(gainedCNEs, lostCNEs)
colnames(total) = c("chr","start","end","id")
colnames(transposons) = c("chr","start","end","id")
#simpleRepeats = rbind(simpleRepeats, transposons)
colnames(simpleRepeats) = c("chr","start","end","id")

goodChrOrder = paste('chr', c(1:22,'X','Y'), sep='')
total$chr = factor(total$chr, levels = goodChrOrder)
total$id = as.factor(total$id)

length(transposons$start)
length(total$start)

#######
setwd('/Users/osipova/Documents/LabDocs/RepeatFiller_docs')
gainedCNEs = read.table('len_gainedCNEs_30bp.bed')[,1:3]
allCNEs = read.table('cons.noExon_noGF_30bp.bed')[,1:3]
transposons = read.table('allTransposons.hg38.bed', sep='\t', header=F)
colnames(transposons) = c("chr","start","end", "id")

total = rbind(data.frame(type = 'G', value = gainedCNEs),
              data.frame(type = 'L', value = allCNEs))
colnames(total) = c("id", "chr","start","end")
goodChrOrder = paste('chr', c(1:22,'X','Y'), sep='')
total$chr = factor(total$chr, levels = goodChrOrder)
total$id = as.factor(total$id)

###### Plot distributions along hg38 chromosomes
library(ggplot2)
#  facet_wrap(~ chr,ncol=2) + # seperate plots for each chr, x-scales can differ from chr to chr

for (i in goodChrOrder){
  print(i)
  chri = subset(total, total$chr==i)
  transi = subset(transposons, transposons$chr == i)
  #simplei = subset(simpleRepeats, simpleRepeats$chr == i)
  cneTotal = ggplot(chri) +
    # histograms
    #geom_histogram(data=subset(chri,chri$id == 'G'), aes(x=subset(chri$start,chri$id == 'G')), aes(y=..density..), fill='red',alpha=0.3,binwidth=1e6) + 
    #geom_histogram(data=subset(chri,chri$id == 'L'), aes(x=subset(chri$start,chri$id == 'L')), aes(y=..density..), fill='blue',alpha=0.3,binwidth=1e6) +
    #geom_histogram(data=transi, aes(x=transi$start),aes(y=..density..), fill='green',alpha=0.3,binwidth=1e6) + 
    
    # density plots
    geom_density(data=subset(chri,chri$id == 'G'), aes(x=subset(chri$start,chri$id == 'G')), adjust=0.3, fill='red',alpha=0.3) +
    geom_density(data=subset(chri,chri$id == 'L'), aes(x=subset(chri$start,chri$id == 'L')), adjust=0.3, fill='blue',alpha=0.43) +
    geom_density(data=transi, aes(x=transi$start), adjust=0.3, fill='green',alpha=0.3) +
    #geom_density(data=simplei, aes(x=simplei$start), adjust=0.3, fill='green',alpha=0.3) +
    
    xlab("Position in the genome") + 
    ylab("count") +
    ggtitle("Density of gained and all CNEs across hg38 chroms")
  png(paste(i, "_gained_vs_all_density_1e6.png", sep=''),1800,500)
  print(cneTotal)
  dev.off()
}
