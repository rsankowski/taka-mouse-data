library(tidyverse)
load("data/tsne_data2.Rda")
source('~/Documents/Single cell analysis/Advanced-plots/20181025-sankowski-et-al-functions.R')

gene_means <- tsne_data2[ ,c(11,14:ncol(tsne_data2))] %>%
  group_by(RaceID_cluster) %>%
  summarise_all(funs(mean))

genes2 <- as.matrix(gene_means[-1])
rownames(genes2) <- gene_means$RaceID_cluster
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

for (i in plot_genes) {
  tryCatch({
    svg(paste0('plots/others/', i, '-violin-plot-linear.svg'), width = 8.57, height = 5.79)
    pl <- violin_plot2(gene=c(i), point_size = 3, logsc = T, .colors = c(colors_many,colors_pat))
    print(pl)
    dev.off()
    
    ggsave(paste0('plots/others/', i, '-violin-plot-logsc.pdf'), width = 8.57, height = 5.79)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  on.exit(dev.off())
} 
