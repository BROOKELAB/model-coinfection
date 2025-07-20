library(tidyverse)

load("runs/rdata/counts_list.RData")
load("runs/rdata/fitness_list.RData")

#calc max fitness
all.counts0 <- map(counts.list, function(x) x$mu0)
all.counts01 <- map(counts.list, function(x) x$mu01)
all.counts05 <- map(counts.list, function(x) x$mu05)
all.counts09 <- map(counts.list, function(x) x$mu09)
all.counts1 <- map(counts.list, function(x) x$mu1)

calc.max.fitness <- function(counts0,countsx,fit){
  fit.0 <- fit
  max.list.0 <- list()
  for(i in seq_along(fit.0$A1B1C1)){
    fit.0[i,][which(counts0[i,] == 0)] <- 0
    max.list.0[[i]] <- max(fit.0[i,])
  }
  max.list.0 <- unlist(max.list.0)
  
  fit.x <- fit
  max.list.x <- list()
  for(i in seq_along(fit.x$A1B1C1)){
    fit.x[i,][which(countsx[i,] == 0)] <- 0
    max.list.x[[i]] <- max(fit.x[i,])
  }
  max.list.x <- unlist(max.list.x)
  
  max.list <- list("mu0" = max.list.0, "mu_plus" = max.list.x)
  return(max.list)
}

apply.max.fitness <- function(counts0, countsx, fit){
  fit.list <- list("noepi" = NA, "epi"= NA, "midepi" = NA,"pos"=NA,"neg" = NA)
  for(i in seq_along(counts0)){
    fit.list[[i]] <- calc.max.fitness(counts0[[i]], countsx[[i]],fit[[i]])
  }
  return(fit.list)
}

max.fitness.01 <- apply.max.fitness(all.counts0, all.counts01, fitness.list)
max.fitness.05 <- apply.max.fitness(all.counts0, all.counts05, fitness.list)
max.fitness.09 <- apply.max.fitness(all.counts0, all.counts09, fitness.list)
max.fitness.1 <- apply.max.fitness(all.counts0, all.counts1, fitness.list)

max.fitness.list <- list("mu01" = max.fitness.01, "mu05" = max.fitness.05,
                         "mu09" = max.fitness.09, "mu1" = max.fitness.1)

#plot max fitness#

#no epistasis
noepi <- map(max.fitness.list, function(x) bind_rows(x$noepi))
noepi.gather <- map(noepi, gather, value = "fitness", key = "mu")

noepi.gather$mu01$mu <- factor(noepi.gather$mu01$mu,levels = c("mu0","mu_plus"))
ggplot(noepi.gather$mu01, aes(x = fitness, fill = mu))+
  geom_density(alpha = 0.5)+
  xlab("Maximum fitness")+
  ylab("Density")+
  scale_x_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1.0,1.5))+
  scale_fill_manual(values = c("grey", "#FED976"))+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 19))+
  guides(fill=guide_legend(title="\U003BC"))
ggsave("runs/figs/max_fitness/noepi_mu01.png", width = 6, height = 4, dpi = 300)

noepi.gather$mu05$mu <- factor(noepi.gather$mu05$mu,levels = c("mu0","mu_plus"))
ggplot(noepi.gather$mu05, aes(x = fitness, fill = mu))+
  geom_density(alpha = 0.5)+
  xlab("Maximum fitness")+
  ylab("Density")+
  scale_x_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1.0,1.5))+
  scale_fill_manual(values = c("grey","#FD8D3C"))+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 19))+
  guides(fill=guide_legend(title="\U003BC"))
ggsave("runs/figs/max_fitness/noepi_mu05.png", width = 6, height = 4, dpi = 300)

noepi.gather$mu09$mu <- factor(noepi.gather$mu09$mu,levels = c("mu0","mu_plus"))
ggplot(noepi.gather$mu09, aes(x = fitness, fill = mu))+
  geom_density(alpha = 0.5)+
  xlab("Maximum fitness")+
  ylab("Density")+
  scale_x_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1.0,1.5))+
  scale_fill_manual(values = c("grey","#E31A1C"))+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 19))+
  guides(fill=guide_legend(title="\U003BC"))
ggsave("runs/figs/max_fitness/noepi_mu09.png", width = 6, height = 4, dpi = 300)

noepi.gather$mu1$mu <- factor(noepi.gather$mu1$mu,levels = c("mu0","mu_plus"))
ggplot(noepi.gather$mu1, aes(x = fitness, fill = mu))+
  geom_density(alpha = 0.5)+
  xlab("Maximum fitness")+
  ylab("Density")+
  scale_x_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1.0,1.5))+
  scale_fill_manual(values = c("grey","#800026"))+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 19))+
  guides(fill=guide_legend(title="\U003BC"))
ggsave("runs/figs/max_fitness/noepi_mu1.png", width = 6, height = 4, dpi = 300)

wilcox.test(noepi$mu01$mu0, noepi$mu01$mu_plus, paired =F) #0.03603
wilcox.test(noepi$mu05$mu0, noepi$mu05$mu_plus, paired =F) #0.0009834
wilcox.test(noepi$mu09$mu0, noepi$mu09$mu_plus, paired =F) #0.0005055
wilcox.test(noepi$mu1$mu0, noepi$mu1$mu_plus, paired =F) #0.001025


#epistasis
epi <- map(max.fitness.list, function(x) bind_rows(x$epi))
epi.gather <- map(epi, gather, value = "fitness", key = "mu")

epi.gather$mu01$mu <- factor(epi.gather$mu01$mu,levels = c("mu0","mu_plus"))
ggplot(epi.gather$mu01, aes(x = fitness, fill = mu))+
  geom_density(alpha = 0.5)+
  xlab("Maximum fitness")+
  ylab("Density")+
  scale_x_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1.0,1.5))+
  scale_fill_manual(values = c("grey", "#FED976"))+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 19))+
  guides(fill=guide_legend(title="\U003BC"))
ggsave("runs/figs/max_fitness/epi_mu01.png",width = 6, height = 4, dpi = 300)

epi.gather$mu05$mu <- factor(epi.gather$mu05$mu,levels = c("mu0","mu_plus"))
ggplot(epi.gather$mu05, aes(x = fitness, fill = mu))+
  geom_density(alpha = 0.5)+
  xlab("Maximum fitness")+
  ylab("Density")+
  scale_x_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1.0,1.5))+
  scale_fill_manual(values = c("grey","#FD8D3C"))+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 19))+
  guides(fill=guide_legend(title="\U003BC"))
ggsave("runs/figs/max_fitness/epi_mu05.png",width = 6, height = 4, dpi = 300)

epi.gather$mu09$mu <- factor(epi.gather$mu09$mu,levels = c("mu0","mu_plus"))
ggplot(epi.gather$mu09, aes(x = fitness, fill = mu))+
  geom_density(alpha = 0.5)+
  xlab("Maximum fitness")+
  ylab("Density")+
  scale_x_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1.0,1.5))+
  scale_fill_manual(values = c("grey","#E31A1C"))+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 19))+
  guides(fill=guide_legend(title="\U003BC"))
ggsave("runs/figs/max_fitness/epi_mu09.png",width = 6, height = 4, dpi = 300)

epi.gather$mu1$mu <- factor(epi.gather$mu1$mu,levels = c("mu0","mu_plus"))
ggplot(epi.gather$mu1, aes(x = fitness, fill = mu))+
  geom_density(alpha = 0.5)+
  xlab("Maximum fitness")+
  ylab("Density")+
  scale_x_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1.0,1.5))+
  scale_fill_manual(values = c("grey","#800026"))+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 19))+
  guides(fill=guide_legend(title="\U003BC"))
ggsave("runs/figs/max_fitness/epi_mu1.png",width = 6, height = 4, dpi = 300)

wilcox.test(epi$mu01$mu0, epi$mu01$mu_plus) #0.239
wilcox.test(epi$mu05$mu0, epi$mu05$mu_plus) #0.01481
wilcox.test(epi$mu09$mu0, epi$mu09$mu_plus) #0.02806
wilcox.test(epi$mu1$mu0, epi$mu1$mu_plus) #0.01801



