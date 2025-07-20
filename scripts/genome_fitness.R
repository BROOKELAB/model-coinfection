library(tidyverse)
library(ggpubr)
library(dunn.test)
library(here)

#calculate fitness for each genotype
load("rdata/counts_list.RData")
load("rdata/gen_lists.RData")
load("rdata/combo.RData")

calc.fitness <- function(gen){
  fitness <- as.data.frame(matrix(ncol = 8, nrow = 1))
  colnames(fitness) <- combo
  fitness[1] <- min(c((gen$i1[1,1] * (1 - abs(gen$i1[1,2] - gen$i1[2,2]))),
                      (gen$i1[2,1] * (1 - abs(gen$i1[2,2] - gen$i1[1,2]))), 
                      (gen$i1[3,1]))) #A1B1C1
  fitness[2] <- min(c((gen$i1[1,1] * (1 - abs(gen$i1[1,2] - gen$i1[2,2]))),
                      (gen$i1[2,1] * (1 - abs(gen$i1[2,2] - gen$i1[1,2]))), 
                      (gen$i2[3,1]))) #A1B1C2
  fitness[3] <- min(c((gen$i1[1,1] * (1 - abs(gen$i1[1,2] - gen$i2[2,2]))),
                      (gen$i2[2,1] * (1 - abs(gen$i2[2,2] - gen$i1[1,2]))), 
                      (gen$i1[3,1]))) #A1B2C1
  fitness[4] <- min(c((gen$i1[1,1] * (1 - abs(gen$i1[1,2] - gen$i2[2,2]))),
                      (gen$i2[2,1] * (1 - abs(gen$i2[2,2] - gen$i1[1,2]))), 
                      (gen$i2[3,1]))) #A1B2C2
  fitness[5] <- min(c((gen$i2[1,1] * (1 - abs(gen$i2[1,2] - gen$i1[2,2]))),
                      (gen$i1[2,1] * (1 - abs(gen$i1[2,2] - gen$i2[1,2]))), 
                      (gen$i1[3,1]))) #A2B1C1
  fitness[6] <- min(c((gen$i2[1,1] * (1 - abs(gen$i2[1,2] - gen$i1[2,2]))),
                      (gen$i1[2,1] * (1 - abs(gen$i1[2,2] - gen$i2[1,2]))), 
                      (gen$i2[3,1]))) #A2B1C2
  fitness[7] <- min(c((gen$i2[1,1] * (1 - abs(gen$i2[1,2] - gen$i2[2,2]))),
                      (gen$i2[2,1] * (1 - abs(gen$i2[2,2] - gen$i2[1,2]))), 
                      (gen$i1[3,1]))) #A2B2C1
  fitness[8] <- min(c((gen$i2[1,1] * (1 - abs(gen$i2[1,2] - gen$i2[2,2]))),
                      (gen$i2[2,1] * (1 - abs(gen$i2[2,2] - gen$i2[1,2]))), 
                      (gen$i2[3,1]))) #A2B2C2
  return(fitness)
}

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

genotype.noepi.fitness <- map_dfr(gen.lists, calc.fitness.noepi)
genotype.epi.fitness <- map_dfr(gen.lists, calc.fitness)
genotype.midepi.fitness <- map_dfr(gen.lists, calc.fitness.midepi)

fitness.list <- list("noepi" = genotype.noepi.fitness, "epi" = genotype.epi.fitness,
                     "pos" = genotype.noepi.fitness,"neg" = genotype.noepi.fitness,
                     "noepi_avg" = genotype.noepi.fitness)

save(fitness.list, file = "rdata/fitness_list.RData")


#calculate weighted fitness per population
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
    weighted.list[[i]] <- calc.weighted.fitness(fitnesses, counts[[i]])
  }
  return(weighted.list)
}

weighted.fitness.list <- map2(fitness.list, counts.list, apply.weighted.fitness)

#mu1/0.9/0.5/0.1 vs mu0 fitness
#without epistasis
ratios.noepi <- bind_cols(weighted.fitness.list$noepi)%>%
  mutate("vs_mu01" = mu01/mu0)%>%
  mutate("vs_mu05" = mu05/mu0)%>%
  mutate("vs_mu09" = mu09/mu0)%>%
  mutate("vs_mu1" = mu1/mu0)%>%
  select(vs_mu01, vs_mu05, vs_mu09, vs_mu1)

ratios.noepi.gather <- gather(ratios.noepi, key = "comparison", value = "ratio")

ggplot(data = ratios.noepi.gather[which(ratios.noepi.gather$comparison == "vs_mu01"),], 
       aes(x = ratio))+
  geom_histogram(alpha = 0.8, fill = "#FED976", color = "white", bins = 40)+
  scale_x_log10(limits = c(0.05, 100000),oob = scales::oob_keep)+
  xlab("Fitness ratio")+
  ylim(c(0,80))+
  ylab("Count")+
  geom_vline(xintercept = 1, linetype = 2)+
  ggtitle("\U003BC = 0.1 vs \U003BC = 0.0")+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 22))
ggsave("figs/fitness_comps/histograms/mu01_vs_mu0_noepi.png",
       width = 5, height = 4, dpi = 300)

ggplot(data = ratios.noepi.gather[which(ratios.noepi.gather$comparison == "vs_mu05"),], 
       aes(x = ratio))+
  geom_histogram(alpha = 0.7, fill = "#FD8D3C", color = "white", bins = 40)+
  scale_x_log10(limits = c(0.05, 100000), oob = scales::oob_keep)+
  xlab("Fitness ratio")+
  ylim(c(0,80))+
  ylab("Count")+
  geom_vline(xintercept = 1, linetype = 2)+
  ggtitle("\U003BC = 0.5 vs \U003BC = 0.0")+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 22))
ggsave("figs/fitness_comps/histograms/mu05_vs_mu0_noepi.png",
       width = 5, height = 4, dpi = 300)

ggplot(data = ratios.noepi.gather[which(ratios.noepi.gather$comparison == "vs_mu09"),], 
       aes(x = ratio))+
  geom_histogram(alpha = 0.6, fill = "#E31A1C", color = "white", bins = 40)+
  scale_x_log10(limits = c(0.05, 100000), oob = scales::oob_keep)+
  xlab("Fitness ratio")+
  ylim(c(0,80))+
  ylab("Count")+
  geom_vline(xintercept = 1, linetype = 2)+
  ggtitle("\U003BC = 0.9 vs \U003BC = 0.0")+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 22))
ggsave("figs/fitness_comps/histograms/mu09_vs_mu0_noepi.png",
       width = 5, height = 4, dpi = 300)

ggplot(data = ratios.noepi.gather[which(ratios.noepi.gather$comparison == "vs_mu1"),], 
       aes(x = ratio))+
  geom_histogram(alpha = 0.5, fill = "#800026", color = "white", bins = 40)+
  scale_x_log10(limits = c(0.05, 100000), oob = scales::oob_keep)+
  xlab("Fitness ratio")+
  ylim(c(0,80))+
  ylab("Count")+
  geom_vline(xintercept = 1, linetype = 2)+
  ggtitle("\U003BC = 1.0 vs \U003BC = 0.0")+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 22))
ggsave("figs/fitness_comps/histograms/mu1_vs_mu0_noepi.png",
       width = 5, height = 4, dpi = 300)

#with epistasis
ratios.epi <- bind_cols(weighted.fitness.list$epi) %>%
  mutate("vs_mu01" = mu01/mu0)%>%
  mutate("vs_mu05" = mu05/mu0)%>%
  mutate("vs_mu09" = mu09/mu0)%>%
  mutate("vs_mu1" = mu1/mu0)%>%
  select(vs_mu01, vs_mu05, vs_mu09, vs_mu1)

ratios.epi.gather <- gather(ratios.epi, key = "comparison", value = "ratio")

ggplot(data = ratios.epi.gather[which(ratios.epi.gather$comparison == "vs_mu01"),], 
       aes(x = ratio))+
  geom_histogram(alpha = 0.8, fill = "#FED976", color = "white", bins = 40)+
  scale_x_log10(limits = c(0.05, 100000), oob = scales::oob_keep)+
  xlab("Fitness ratio")+
  ylim(c(0,60))+
  ylab("Count")+
  geom_vline(xintercept = 1, linetype = 2)+
  ggtitle("\U003BC = 0.1 vs \U003BC = 0.0")+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 22))
ggsave("figs/fitness_comps/histograms/mu01_vs_mu0_epi.png",
       width = 5, height = 4, dpi = 300)

ggplot(data = ratios.epi.gather[which(ratios.epi.gather$comparison == "vs_mu05"),], 
       aes(x = ratio))+
  geom_histogram(alpha = 0.7, fill = "#FD8D3C", color = "white", bins = 40)+
  scale_x_log10(limits = c(0.05, 100000), oob = scales::oob_keep)+
  xlab("Fitness ratio")+
  ylim(c(0,60))+
  ylab("Count")+
  geom_vline(xintercept = 1, linetype = 2)+
  ggtitle("\U003BC = 0.5 vs \U003BC = 0.0")+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 22))
ggsave("figs/fitness_comps/histograms/mu05_vs_mu0_epi.png",
       width = 5, height = 4, dpi = 300)

ggplot(data = ratios.epi.gather[which(ratios.epi.gather$comparison == "vs_mu09"),], 
       aes(x = ratio))+
  geom_histogram(alpha = 0.6, fill = "#E31A1C", color = "white", bins = 40)+
  scale_x_log10(limits = c(0.05, 100000), oob = scales::oob_keep)+
  xlab("Fitness ratio")+
  ylim(c(0,60))+
  ylab("Count")+
  geom_vline(xintercept = 1, linetype = 2)+
  ggtitle("\U003BC = 0.9 vs \U003BC = 0.0")+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 22))
ggsave("figs/fitness_comps/histograms/mu09_vs_mu0_epi.png",
       width = 5, height = 4, dpi = 300)

ggplot(data = ratios.epi.gather[which(ratios.epi.gather$comparison == "vs_mu1"),], 
       aes(x = ratio))+
  geom_histogram(alpha = 0.5, fill = "#800026", color = "white", bins = 40)+
  scale_x_log10(limits = c(0.05, 100000), oob = scales::oob_keep)+
  xlab("Fitness ratio")+
  ylim(c(0,60))+
  ylab("Count")+
  geom_vline(xintercept = 1, linetype = 2)+
  ggtitle("\U003BC = 1.0 vs \U003BC = 0.0")+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 22))
ggsave("figs/fitness_comps/histograms/mu1_vs_mu0_epi.png",
       width = 5, height = 4, dpi = 300)


#positive dominance
ratios.pos <- bind_cols(weighted.fitness.list$pos) %>%
  mutate("vs_mu01" = mu01/mu0)%>%
  mutate("vs_mu05" = mu05/mu0)%>%
  mutate("vs_mu09" = mu09/mu0)%>%
  mutate("vs_mu1" = mu1/mu0)%>%
  select(vs_mu01, vs_mu05, vs_mu09, vs_mu1)

ratios.pos.gather <- gather(ratios.pos, key = "comparison", value = "ratio")

ggplot(data = ratios.pos.gather[which(ratios.pos.gather$comparison == "vs_mu1"),], 
       aes(x = ratio))+
  geom_histogram(alpha = 0.5, fill = "#800026", color = "white", bins = 40)+
  scale_x_log10(limits = c(0.05, 100000), oob = scales::oob_keep)+
  xlab("Fitness ratio")+
  ylim(c(0,40))+
  ylab("Count")+
  geom_vline(xintercept = 1, linetype = 2)+
  ggtitle("\U003BC = 1.0 vs \U003BC = 0.0")+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 22))
ggsave("figs/fitness_comps/histograms/mu1_vs_mu0_pos.png",
       width = 5, height = 4, dpi = 300)

#negative dominance
ratios.neg <- bind_cols(weighted.fitness.list$neg) %>%
  mutate("vs_mu01" = mu01/mu0)%>%
  mutate("vs_mu05" = mu05/mu0)%>%
  mutate("vs_mu09" = mu09/mu0)%>%
  mutate("vs_mu1" = mu1/mu0)%>%
  select(vs_mu01, vs_mu05, vs_mu09, vs_mu1)

ratios.neg.gather <- gather(ratios.neg, key = "comparison", value = "ratio")

ggplot(data = ratios.neg.gather[which(ratios.neg.gather$comparison == "vs_mu1"),], 
       aes(x = ratio))+
  geom_histogram(alpha = 0.5, fill = "#800026", color = "white", bins = 40)+
  scale_x_log10(limits = c(0.05, 100000), oob = scales::oob_keep)+
  xlab("Fitness ratio")+
  ylim(c(0,40))+
  ylab("Count")+
  geom_vline(xintercept = 1, linetype = 2)+
  ggtitle("\U003BC = 1.0 vs \U003BC = 0.0")+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 22))
ggsave("figs/fitness_comps/histograms/mu1_vs_mu0_neg.png",
       width = 5, height = 4, dpi = 300)

#cross-segment avg
ratios.noepi.avg <- bind_cols(weighted.fitness.list$noepi_avg)%>%
  mutate("vs_mu01" = mu01/mu0)%>%
  mutate("vs_mu05" = mu05/mu0)%>%
  mutate("vs_mu09" = mu09/mu0)%>%
  mutate("vs_mu1" = mu1/mu0)%>%
  select(vs_mu01, vs_mu05, vs_mu09, vs_mu1)

ratios.noepi.avg.gather <- gather(ratios.noepi.avg, key = "comparison", value = "ratio")

ggplot(data = ratios.noepi.avg.gather[which(ratios.noepi.avg.gather$comparison == "vs_mu1"),], 
       aes(x = ratio))+
  geom_histogram(alpha = 0.5, fill = "#800026", color = "white", bins = 40)+
  scale_x_log10(limits = c(0.05, 100000), oob = scales::oob_keep)+
  xlab("Fitness ratio")+
  ylim(c(0,40))+
  ylab("Count")+
  geom_vline(xintercept = 1, linetype = 2)+
  ggtitle("\U003BC = 1.0 vs \U003BC = 0.0")+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 22))
ggsave("figs/fitness_comps/histograms/mu1_vs_mu0_noepi_avg.png",
       width = 5, height = 4, dpi = 300)

