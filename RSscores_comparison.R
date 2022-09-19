setwd('/Users/osipova/Documents/R_scripts/May_2019_CNE_analysis')

gapFill_means = scan('all_chrom_gapFill_RS_means.lst')
noGapFill_means = scan('all_chrom_noGapFill_RS_means.lst')

abs_RS_scores_gF_nGF = rbind(data.frame(type = 'noGapFill',value = noGapFill_means),
          data.frame(type = 'gapFill',value = gapFill_means))

library(ggplot2)
p1 = ggplot(abs_RS_scores_gF_nGF, aes(x = type, y = value, fill = type)) + geom_violin() + geom_boxplot(width=0.1)
p2 = p1 + labs(y = 'RS score GERP', title = 'abs GERP RS score per bp', x = '')
p3 = p2 + scale_fill_manual(values=c('gapFill' = "coral", 'noGapFill' = "gray45")) + theme_classic()
p3 + stat_summary(fun.y=median, geom="point", size=1, color="black")

RS_combined_gF = read.table('all_chroms_combined_RS_gained.lst')
RS_combined_nGF = read.table('all_chroms_combined_RS_lost.lst')
delta_RS_gF = RS_combined_gF$V1 - RS_combined_gF$V2
delta_RS_nGF = RS_combined_nGF$V2 - RS_combined_nGF$V1

delta_RS_scores_gF_nGF = rbind(data.frame(type = 'lost', value = delta_RS_nGF),
                               data.frame(type = 'gained', value = delta_RS_gF))

p1 = ggplot(delta_RS_scores_gF_nGF, aes(x = type, y = value, fill = type)) + geom_violin() 
p2 = p1 + geom_boxplot(width=0.1, outlier.shape = NA) + labs(y = 'difference in rejected substitutions scores per base pair', x = '')
p3 = p2 + scale_fill_manual(values=c('gained' = "coral", 'lost' = "gray45")) + theme_classic()
#p3 + stat_summary(fun.y=median, geom="point", size=1, color="black")
p3 + scale_y_continuous(limits = c(-6, 8)) + coord_flip()


### exon lengths distribution ###
exonLenshg38 = scan('exonLen.hg38.ensGene_exon.bed')
cdsLenshg38 = scan('cdsLen.hg38.ensGene_cds.bed')
exonLensmm10 = scan('exonLen.mm10.ensGene_exon.bed')
cdsLensmm10 = scan('cdsLen.mm10.ensGene_cds.bed')

par(mfrow=c(2,2))
hist(head(sort(exonLenshg38),-70000), breakes=5000)
hist(head(sort(cdsLenshg38),-10000), breakes=5000)
hist(head(sort(exonLensmm10),-50000), breakes=5000)
hist(head(sort(cdsLensmm10),-10000), breakes=5000)
