library(tidyverse)

#viable genomes stats (immune resistance)
load("rdata/resistance_viables_list.RData")

viables.15.ab <- map_dfc(viables.15.noepi, function(x) x$ab)
viables.15.ab <- gather(viables.15.ab, key = "mu", value = "viables")
viables.15.c <- map_dfc(viables.15.noepi, function(x) x$c)
viables.15.c <- gather(viables.15.c, key = "mu", value = "viables")

mu.cols <- c("black", RColorBrewer::brewer.pal(9,"YlOrRd")[c(3,5,7,9)])
ggplot(viables.15.ab, aes(x = mu, y = viables, fill = mu))+
  geom_violin()+
  ylab("Viable genomes")+
  scale_y_log10(limits = c(10,1000))+
  xlab("\U003BC")+
  scale_x_discrete(labels=c("mu0" = "0.0", "mu01" = "0.1","mu05" = "0.5", 
                            "mu09" = "0.9","mu1" = "1.0"))+
  scale_fill_manual(values = mu.cols)+
  ggtitle("AB < C")+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 22),
        legend.position = "none")
ggsave("figs/resistance/noepi_15_ab.png",
       width = 5, height = 4)


ggplot(viables.15.c, aes(x = mu, y = viables, fill = mu))+
  geom_violin()+
  ylab("Viable genomes")+
  scale_y_log10(limits = c(10,1000))+
  xlab("\U003BC")+
  scale_x_discrete(labels=c("mu0" = "0.0", "mu01" = "0.1","mu05" = "0.5", 
                            "mu09" = "0.9","mu1" = "1.0"))+
  scale_fill_manual(values = mu.cols)+
  ggtitle("C < AB")+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 22),
        legend.position = "none")
ggsave("figs/resistance/noepi_15_c.png",
       width = 5, height = 4)

#calculate correlations
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

calc.viable.cor(viables.5.ab) #r = na #p = na
calc.viable.cor(viables.5.c) # r = -0.03925779 #p = 0.3810456

calc.viable.cor(viables.15.ab) #r = 0.08414793 #p = 0.06007655
calc.viable.cor(viables.15.c) # r = 0.05059689 #p = 0.2587828

calc.viable.cor(viables.25.ab) #r = -0.07912437 #p = 0.07712468
calc.viable.cor(viables.25.c) # r = -0.06417774 #p = 0.1518739

calc.viable.cor(viables.35.ab) #r = -0.003346311 #p = 0.9405022
calc.viable.cor(viables.35.c) # r = -0.07171804 #p = 0.109219


#calculate fitness for each genotype
load("rdata/resistance_fitness_list.RData")
load("rdata/combo.RData")

calc.fitness.noepi <- function(gen){
  fitness <- as.data.frame(matrix(ncol = 8, nrow = 1))
  colnames(fitness) <- combo
  fitness[1] <- min(c(gen$i1[1,1],
                      gen$i1[2,1], 
                      gen$i1[3,1])) #A1B1C1
  fitness[2] <- min(c(gen$i1[1,1],
                      gen$i1[2,1], 
                      gen$i2[3,1])) #A1B1C2
  fitness[3] <- min(c(gen$i1[1,1],
                      gen$i2[2,1], 
                      gen$i1[3,1])) #A1B2C1
  fitness[4] <- min(c(gen$i1[1,1],
                      gen$i2[2,1], 
                      gen$i2[3,1])) #A1B2C2
  fitness[5] <- min(c(gen$i2[1,1],
                      gen$i1[2,1], 
                      gen$i1[3,1])) #A2B1C1
  fitness[6] <- min(c(gen$i2[1,1],
                      gen$i1[2,1], 
                      gen$i2[3,1])) #A2B1C2
  fitness[7] <- min(c(gen$i2[1,1],
                      gen$i2[2,1], 
                      gen$i1[3,1])) #A2B2C1
  fitness[8] <- min(c(gen$i2[1,1],
                      gen$i2[2,1], 
                      gen$i2[3,1])) #A2B2C2
  return(fitness)
}
apply.fitness <- function(gen){
  genotype.fitness <- list("abc" = NA, "ab" = NA, "c" = NA)
  for(i in seq_along(gen)){
    genotype.fitness[[i]] <- map_dfr(gen[[i]], calc.fitness.noepi)
  }
  return(genotype.fitness)
}

res.fitness.list <- map(res.fitness.list, apply.fitness)

#get counts for each genotype
load("rdata/resistance_counts_list.RData")

#calculate weighted fitnesses per population
calc.weighted.fitness <- function(fit, counts){
  weighted.fit <- list()
  for(i in seq_along(rownames(fit))){
    weighted.fit[[i]] <- weighted.mean(fit[i,], counts[i,])
  }
  weighted.fit <- unlist(weighted.fit)
  return(weighted.fit)
}

apply.weighted.fitness <- function(fitnesses, counts){
  weighted.list <- list("mu0" = NA, "mu01" = NA, "mu05" = NA, "mu09" = NA, "mu1" = NA)
  for(i in seq_along(counts)){
    weighted.list[[i]] <- map2(fitnesses,counts[[i]], calc.weighted.fitness)
  }
  return(weighted.list)
}

res.weighted.fitness.list <- map2(res.fitness.list, res.counts.list, apply.weighted.fitness)


#per run mu1/mu0 fitness
comp.fitness <- function(weighted.fit){
  ratio.abc <- (weighted.fit$mu1$abc) / (weighted.fit$mu0$abc)
  ratio.ab <- (weighted.fit$mu1$ab) / (weighted.fit$mu0$ab)
  ratio.c <- (weighted.fit$mu1$c) / (weighted.fit$mu0$c)
  ratios <- as.data.frame(cbind("A = B = C" = ratio.abc, "AB < C" = ratio.ab, 
                                "C < AB" = ratio.c))
  return(ratios)
}

res.fitness.comps <- map(res.weighted.fitness.list, comp.fitness)

fitcomp <- res.fitness.comps
names(fitcomp) <- c("5%","15%","25%","35%")
fitcomp <- bind_rows(fitcomp, .id = "percent")
fitcomp.gather <- gather(fitcomp, value = "Fitness ratio", key = "Antigenic scheme", - percent)

scheme.cols <- c("#DEEBF7","#3182BD","#9ECAE1")
fitcomp.gather$percent <- factor(fitcomp.gather$percent, levels = c("5%","15%","25%","35%"))
ggplot(fitcomp.gather, aes(x = percent, y = `Fitness ratio`), color = "black")+
  geom_point(cex = 4, shape = 21, position = position_dodge(width = 0.5),
             aes(fill = `Antigenic scheme`), alpha = 0.3)+
  geom_hline(yintercept = 1, linetype = 2)+
  xlab("Percent resistance")+
  ylab("\U003BC 1.0 vs \U003BC 0.0 fitness ratio")+
  scale_y_log10()+
  scale_fill_manual(values = scheme.cols)+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 19),
        legend.title = element_text(size = 22))+
  guides(fill = guide_legend(override.aes = list(alpha = 1)))
ggsave("figs/resistance/resistance_fitness.png", 
       width = 13, height = 5, dpi = 300)









