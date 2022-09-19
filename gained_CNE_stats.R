setwd('/Users/osipova/Documents/LabDocs/RepeatFiller_docs')
library(ggplot2)

gainedCNEs = read.table('len_gainedCNEs_30bp.bed')
allCNEs = read.table('cons.noExon_noGF_30bp.bed')

colnames(allCNEs) = c('chr', 'start', 'end', 'size')
colnames(gainedCNEs) = colnames(allCNEs)

##### size distribution #####
# combine data keeping the type
combined_CNEs = rbind(data.frame(type = 'gained', value = gainedCNEs$size),
                      data.frame(type = 'noRepeatFill', value = allCNEs$size))

# density plots
ggplot(combined_CNEs, aes(x=value, color=type)) + geom_density(alpha=0.3) + theme_bw() +
  labs(x = 'element size, bp') + xlim(30, 200) + scale_color_manual(values=c("#F47E53", "#999999")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#####
# density plots: double log scale
log_gainedCNEs = log(gainedCNEs$size)
log_allCNEs = log(allCNEs$size)

log_combined_CNEs = rbind(data.frame(type = 'gained', value = log_gainedCNEs),
                      data.frame(type = 'noRepeatFill', value = log_allCNEs))

ggplot(log_combined_CNEs, aes(x=value, color=type)) + geom_line(stat="density") + theme_bw() +
  labs(x = 'element size, bp') + scale_y_log10() + xlim(3,8) + 
  scale_color_manual(values=c("#F47E53", "#999999")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())


#####
count_gainedTab = as.data.frame(table(gainedCNEs$size), stringsAsFactors=F)
count_allTab = as.data.frame(table(allCNEs$size), stringsAsFactors=F)

# rename columns
colnames(count_gainedTab) = c('size', 'count')
colnames(count_allTab) = colnames(count_gainedTab)

# fix the type of the size column
count_gainedTab = cbind(as.integer(count_gainedTab$size), count_gainedTab$count)
count_allTab = cbind(as.integer(count_allTab$size), count_allTab$count)

# combine count tables
count_combineTab = rbind(data.frame(type = 'gained', size = count_gainedTab[,1], count = count_gainedTab[,2]),
                         data.frame(type = 'noRepeatFill', size = count_allTab[,1], count = count_allTab[,2]))

ggplot(count_combineTab, aes(x=size, y=count, color=type)) + geom_point(size=0.6) + theme_bw() + 
  scale_color_manual(values=c("#F47E53", "#999999")) + scale_y_log10() + scale_x_log10() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + labs(x='element size, bp', y='element count')
  


#####
# violin plots with boxplots
ggplot(combined_CNEs, aes(x=type, y=value, fill=type)) + geom_violin() + ylim(30,200) + theme_bw() + 
  labs(y = 'element size, bp', x='') + scale_fill_manual(values=c("#F47E53", "#999999")) + 
  geom_boxplot(width=0.05, outlier.shape = NA) + theme(panel.grid.major = element_blank(),
                                                       panel.grid.minor = element_blank()) + coord_flip()

##### the difference in sizes is statistically significant 
wilcox.test(gainedCNEs$size, allCNEs$size)
# W = 1.1886e+10, p-value < 2.2e-16

##### AT percent stats #####
at_gainedCNEs = scan('gainedCNEs_30bp.ATpercent.txt')
at_allCNEs = scan('enumerated_CNEs_noRepeatFill_30bp.ATpercent.txt')
at_shuffled = scan('1.out.shuffled.ATpercent.txt')

# combine data keeping the type
at_combined_CNEs = rbind(data.frame(type = 'gained', value = at_gainedCNEs),
                         data.frame(type = 'noRepeatFill', value = at_allCNEs),
                         data.frame(type = 'shuffled', value = at_shuffled))
# violin plots with boxplots
ggplot(at_combined_CNEs, aes(x=type, y=value, fill=type)) + geom_violin() + theme_bw() + 
  labs(y = 'AT content, %', x='') + scale_fill_manual(values=c("#F47E53", "#999999", "#DDDDDD")) +
  geom_boxplot(width=0.05, outlier.shape = NA) + theme(legend.position = "none",
            panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##### all the differences are significant #####
wilcox.test(at_gainedCNEs, at_allCNEs)
# W = 2.0498e+10, p-value < 2.2e-16
wilcox.test(at_gainedCNEs, at_shuffled)
# W = 2.5807e+10, p-value < 2.2e-16
wilcox.test(at_allCNEs, at_shuffled)
# W = 7.8286e+11, p-value < 2.2e-16

#####
# at_content_gained_vs_all