library(tidyverse)

#load viable genome stats
load("rdata/viables_list.RData")

viables.list <- map(viables.list, bind_cols)
viables.list <- map(viables.list, gather, key = "mu", value = "viables")

#viable genome violin plots
plot.violin <- function(viable){
  mu.cols <- c("black", RColorBrewer::brewer.pal(9,"YlOrRd")[c(3,5,7,9)])
  vio.plot <- ggplot(data = viable, aes(x = mu, y = viables, fill = mu))+
    geom_violin()+
    ylab("Viable genomes")+
    scale_y_log10()+
    xlab("\U003BC")+
    scale_x_discrete(labels=c("mu0" = "0.0", "mu01" = "0.1","mu05" = "0.5", 
                              "mu09" = "0.9","mu1" = "1.0"))+
    scale_fill_manual(values = mu.cols)+
    theme_bw()+
    theme(axis.text = element_text(size = 19),
          axis.title = element_text(size = 22),
          plot.title = element_text(size = 22),
          legend.position = "none")
  
  return(vio.plot)
}

viable.violins <- map(viables.list, plot.violin)

noepi <- viable.violins$noepi
noepi + ggtitle("No epistasis") + theme(plot.title = element_text(size = 22))
ggsave("figs/violins/noepi_viables.png", width = 5, height = 4, dpi = 300)

epi <- viable.violins$epi
epi + ggtitle("Epistasis") + theme(plot.title = element_text(size = 22))
ggsave("figs/violins/epi_viables.png", width = 5, height = 4, dpi = 300)

pos <- viable.violins$pos
pos + ggtitle("Positive dominance") + theme(plot.title = element_text(size = 22))
ggsave("figs/violins/pos_viables.png",width = 5, height = 4, dpi = 300)


neg <- viable.violins$neg
neg + ggtitle("Negative dominance") + theme(plot.title = element_text(size = 22))
ggsave("figs/violins/neg_viables.png",width = 5, height = 4, dpi = 300)

avg <- viable.violins$avg
avg + ggtitle("Cross-segment averaged") + theme(plot.title = element_text(size = 22))
ggsave("figs/violins/avg.png",width = 5, height = 4, dpi = 300)

#calculate mu vs viability correlations
calc.viable.cor <- function(viable){
  cor.mat <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 2))
  colnames(cor.mat) <- c("r","p")
  viable$mu[which(viable$mu == "mu0")] <- 0
  viable$mu[which(viable$mu == "mu01")] <- 0.1
  viable$mu[which(viable$mu == "mu05")] <- 0.5
  viable$mu[which(viable$mu == "mu09")] <- 0.9
  viable$mu[which(viable$mu == "mu1")] <- 1
  viable$mu <- as.numeric(viable$mu)
  cor.mat[1,] <- c((cor.test(viable$mu, viable$viables, method = "pearson"))$estimate,
                   (cor.test(viable$mu, viable$viables, method = "pearson"))$p.value)
  
  return(cor.mat)
}

viable.cors <- map(viables.list, calc.viable.cor)
viable.cor.mat <- bind_rows(viable.cors)
rownames(viable.cor.mat) <- names(viable.cors)

write.csv(viable.cor.mat, "figs/violins/viable_cors.csv")


