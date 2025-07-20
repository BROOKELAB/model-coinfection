library(tidyverse)

load("rdata/counts_list.RData")
load("rdata/fitness_list.RData")

#change in per-genotype gene contribution between parental generation and offspring
freq.check <- function(offspring.count, gen.fit){
  offspring.freq <- as.data.frame(t(offspring.count/(sum(offspring.count))))
  colnames(offspring.freq) <- "freq"
  offspring.tab <- mutate(offspring.freq,
                          "parent1.id" = c(1,0.67,0.67,0.33,0.67,0.33,0.33,0),
                          "parent2.id" = c(0, 0.33,0.33,0.67,0.33,0.67,0.67,1))
  offspring.tab <- mutate(offspring.tab,
                          "freqxid.1" = freq * parent1.id,
                          "freqxid.2" = freq * parent2.id)
  
  offspring.fx <- colSums(offspring.tab)[c(4,5)] - 0.5
  names(offspring.fx) <- c("parent1", "parent2")
  parent.diff <- c("parent1" = unname(gen.fit[1] - gen.fit[8]),
                   "parent2" = unname(gen.fit[8] - gen.fit[1]))
  
  freq.mat <- as.data.frame(bind_rows(offspring.fx, parent.diff))
  rownames(freq.mat) <- c("offspring_fx","parental_diff")
  return(freq.mat)
}

freq.table <- function(offspring.count, gen.fit){
  comp <- list()
  for(i in seq_along(rownames(offspring.count))){
    comp[[i]] <- freq.check(offspring.count[i,], gen.fit[i,])
  }
  comp <- map(comp, function(x) t(x)[1,])
  comp <- bind_rows(comp)
  comp <- mutate(comp,"type" = "mismatch")
  for(i in seq_along(comp$offspring_fx)){
    if((comp$offspring_fx[[i]] < 0) && (comp$parental_diff[[i]] < 0)){
      comp$type[[i]] <- "match"
    }
    if((comp$offspring_fx[[i]] >= 0) && (comp$parental_diff[[i]] >= 0)){
      comp$type[[i]] <- "match"
    }
  }
  return(comp)
}

calc.fx <- function(offspring.count, gen.fit){
  fx <- list("mu0"=NA,"mu01"=NA,"mu05"=NA,"mu09"=NA,"mu1"=NA)
  for(i in seq_along(offspring.count)){
    fx[[i]] <- freq.table(offspring.count[[i]], gen.fit)
  }
  return(fx)
}

fx.list <- map2(counts.list, fitness.list, calc.fx)

mu0.list <- map(fx.list, function(x) x$mu0)
mu1.list <- map(fx.list, function(x) x$mu1)

cor.mu <- function(mu.list){
  cor <- suppressWarnings(cor.test(mu.list$parental_diff, mu.list$offspring_fx, method = "spearman"))
  cor.p <- cor$p.value
  cor.r <- cor$estimate
  cor.results <- as.data.frame(matrix(nrow = 2, ncol = 1))
  rownames(cor.results) <- c("correlation","p_value")
  cor.results[,1] <- c(cor.r, cor.p)
  return(cor.results)
}

cor.mu0 <- map(mu0.list, cor.mu)
cor.mu0 <- map(cor.mu0, setNames, "mu0")

cor.mu1 <- map(mu1.list, cor.mu)
cor.mu1 <- map(cor.mu1, setNames, "mu1")

cor.list <- map2(cor.mu0, cor.mu1, bind_cols)

for(i in seq_along(cor.list)){
  write.csv(cor.list[[i]], file = paste0("runs/figs/generation_comps/stats/",
                                         names(cor.list)[[i]],".csv"))
}

plot.change <- function(mu){
  mu.plots <- ggplot(data = mu, aes(x = parental_diff, y = offspring_fx, color = type))+
    geom_point(cex = 3)+
    xlab("Fitness(genotype i) - fitness(genotype j)")+
    ylab("∆ Population frequency of i")+
    ylim(c(-0.5,0.5))+
    scale_color_manual(values = c("black","red"))+
    theme_bw()+
    theme(axis.title = element_text(size = 22),
          axis.text = element_text(size = 19),
          plot.title = element_text(size = 22),
          legend.position = "none")
  
  return(mu.plots)
}

mu0.plots <- map(mu0.list, plot.change)
mu1.plots <- map(mu1.list, plot.change)

#no epistasis
noepi.mu0 <- mu0.plots$noepi
noepi.mu0 + ggtitle("No epistasis")
ggsave("runs/figs/generation_comps/noepi_mu0.png", width = 7, height = 5, dpi = 300)
noepi.mu1 <- mu1.plots$noepi
noepi.mu1 + ggtitle("No epistasis")
ggsave("runs/figs/generation_comps/noepi_mu1.png", width = 7, height = 5, dpi = 300)

#epistasis
epi.mu0 <- mu0.plots$epi
epi.mu0 + ggtitle("Epistasis")
ggsave("runs/figs/generation_comps/epi_mu0.png", width = 7, height = 5, dpi = 300)
epi.mu1 <- mu1.plots$epi
epi.mu1 + ggtitle("Epistasis")
ggsave("runs/figs/generation_comps/epi_mu1.png", width = 7, height = 5, dpi = 300)

#positive dominance
pos.mu0 <- mu0.plots$pos
pos.mu0 + ggtitle("Positive dominance")
ggsave("runs/figs/generation_comps/pos_mu0.png", width = 7, height = 5, dpi = 300)
pos.mu1 <- mu1.plots$pos
pos.mu1 + ggtitle("Positive dominance")
ggsave("runs/figs/generation_comps/pos_mu1.png", width = 7, height = 5, dpi = 300)

#negative dominance
neg.mu0 <- mu0.plots$neg
neg.mu0 + ggtitle("Negative dominance")
ggsave("runs/figs/generation_comps/neg_mu0.png", width = 7, height = 5, dpi = 300)
neg.mu1 <- mu1.plots$neg
neg.mu1 + ggtitle("Negative dominance")
ggsave("runs/figs/generation_comps/neg_mu1.png", width = 7, height = 5, dpi = 300)



