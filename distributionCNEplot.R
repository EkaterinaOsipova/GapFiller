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


################
# install.packages('VennDiagram')
library(VennDiagram)

grid.newpage()
A = c(1:100371)
B = c(10372:109868)
venn.diagram(list("GapFiller"=A, "noGapFiller"=B), fill = c("pink","light blue"),lty = rep("blank",2),
             alpha = rep(0.5, 2),cex=1.5, filename='Venn.png')

A = c(1:16066)
B = c(48:16147)
venn.diagram(list("GapFiller"=A, "noGapFiller"=B), fill = c("pink","light blue"),lty = rep("blank",2),
             alpha = rep(0.5, 2),cex=1.5, filename='Venn2.png')

# differential losses in hummingbirds A and B
A = c(1:973)
B = c(594:1683)
# Intact genes in A and B hummingbirds
A = c(1:7303)
B = c(2570:8692)
# Intact genes in AB (=HLcalAn2) and C (=calAnn1)
A = c(1:8697)
B = c(1548:9422)
# ABC vs D(G+L+M) comparison
A = c(1:619)
B = c(68:5136)
A = c(1:9277)
B = c(3525:8601)
A = c(1:67)
B = c(3:5748)
venn.diagram(list("A"=A, "B"=B), fill = c("pink","light blue"),lty = rep("blank",2), alpha = rep(0.5, 2),cex=1.5, filename='Venn.png')
#
A = c(1:88)
C = c(14:102)
D = c(10:168)
venn.diagram(list(" "=A, " "=C, " "=D), fill = c("red","blue","green"),lty = rep("blank",3), alpha = rep(0.15, 3),cex=3, filename='Venn.png')
#  
  
grid.newpage()
venn.plot <- draw.single.venn(area = 88, lwd = 5, lty = "blank", cex = 4, fill = "red", alpha = 0.15)
venn.plot <- draw.single.venn(area = 34, lwd = 5, lty = "blank", cex = 4, fill = "blue", alpha = 0.15)
venn.plot <- draw.single.venn(area = 44, lwd = 5, lty = "blank", cex = 4, fill = "green", alpha = 0.15)
venn.plot <- draw.pairwise.venn(area1 = 88,area2 = 34, cross.area = 14,cex=4, fill = c("red","blue"), alpha=0.15)
venn.plot <- draw.triple.venn(area1=88,area2=34,area3=44,n12=14,n23=10,n13=11,n123=9,cex=4, fill=c("red","blue","green"), alpha=0.15,lty="blank")
plot(venn.plot)

png("Venn.png",1800,500)

######
library(VennDiagram)

grid.newpage()
A = c(1:163)
B = c(3:166)
venn.diagram(list(" "=A, "  "=B), fill = c("coral","gray45"),lty = rep("blank",2),
             alpha = rep(0.5, 2),cex=0.1, filename='Venn_gF_nGF.png')


dev.off()
