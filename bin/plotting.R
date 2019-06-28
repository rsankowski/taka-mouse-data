library(tidyverse)
load("data/tsne_data2.Rda")
source('~/Documents/Single cell analysis/Advanced-plots/20181025-sankowski-et-al-functions.R')

tsne_data2 <- tsne_data2[tsne_data2$Age == "16_w",]
tsne_data2$RaceID_cluster <- as.character(tsne_data2$RaceID_cluster)

clust <- as.data.frame(table(tsne_data2$RaceID_cluster))
clust <- clust[clust$Freq>30,]
tsne_data2 <- tsne_data2[tsne_data2$RaceID_cluster %in% clust$Var1,]

gene_means <- tsne_data2[ ,c(11,14:ncol(tsne_data2))] %>%
  group_by(RaceID_cluster) %>%
  summarise_all(funs(mean))

genes2 <- as.matrix(gene_means[1:nrow(gene_means), -1])
rownames(genes2) <- gene_means$RaceID_cluster[1:nrow(gene_means)]
ord_clust <- hclust(dist(genes2))
ord_clust <- rownames(genes2)[ord_clust$order]

#plot
source('~/Documents/Single cell analysis/Advanced-plots/20190618_violin_plot_for_ori_data.R')

plot_genes <- c('Cx3cr1', 'Tgfbr1', 'P2ry12', 'Gpr34', 'Sall1', 'Slc2a5', 'Tmem119', 'Siglech', 'Csf1r', 'Fcrls', 'Hexb')

for (i in plot_genes) {
  tryCatch({
    svg(paste0('plots/others/', i, '-violin-plot-linear.svg'), width = 8.57, height = 5.79)
    pl <- violin_plot2(gene=c(i), point_size = 3, logsc = F, .colors = c(colors_many,colors_pat))
    print(pl)
    dev.off()
    
    ggsave(paste0('plots/others/', i, '-violin-plot-linear.pdf'), width = 8.57, height = 5.79)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  on.exit(dev.off())
} 

#logscale of expressions
for (i in plot_genes) {
  tryCatch({
    svg(paste0('plots/others/', i, '-violin-plot-logsc.svg'), width = 8.57, height = 5.79)
    pl <- violin_plot2(gene=c(i), point_size = 3, .colors = c(colors_many,colors_pat)) + scale_y_log10()
    print(pl)
    dev.off()
    
    ggsave(paste0('plots/others/', i, '-violin-plot-logsc.pdf'), width = 8.57, height = 5.79)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  on.exit(dev.off())
} 

#boxplots
for (i in plot_genes) {
  tryCatch({
    svg(paste0('plots/others/', i, '-box-plot-linear.svg'), width = 8.57, height = 5.79)
    pl <- box_plot2(gene=c(i), point_size = 3, logsc = F, .colors = c(colors_many,colors_pat))
    print(pl)
    dev.off()
    
    ggsave(paste0('plots/others/', i, '-box-plot-linear.pdf'), width = 8.57, height = 5.79)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  on.exit(dev.off())
} 

#logscale of expressions
for (i in plot_genes) {
  tryCatch({
    svg(paste0('plots/others/', i, '-box-plot-logsc.svg'), width = 8.57, height = 5.79)
    pl <- box_plot2(gene=c(i), point_size = 3, .colors = c(colors_many,colors_pat)) + scale_y_log10()
    print(pl)
    dev.off()
    
    ggsave(paste0('plots/others/', i, '-box-plot-logsc.pdf'), width = 8.57, height = 5.79)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  on.exit(dev.off())
} 


tsne_data2$RaceID_cluster <- factor(tsne_data2$RaceID_cluster, levels = ord_clust)
#tsne map of clusters
tsne_plot <- ggplot(na.omit(tsne_data2), aes(V1, V2, fill = RaceID_cluster)) +
  geom_point(pch=21, size=4, stroke=0.25)+
  scale_fill_manual(values=c(colors_many, colors_pat)) +
  theme_void()

tsne_plot

svg(paste0('plots/tsne/clusters-tsne.svg'), width = 8.57, height = 5.79)
tsne_plot
dev.off()

ggsave(paste0('plots/tsne/clusters-tsne.pdf'), width = 8.57, height = 5.79)

