library(tidyverse)

#load viable genome stats
load("runs/next_gen_viables_list.RData")

all.mu0 <- map(ng.viables.list, function(x) unlist(x$mu0))
all.mu0 <- bind_rows(all.mu0, .id = "ng")
all.mu0.gather <- gather(all.mu0, value = "viables", key = "ng")

all.mu1 <- map(ng.viables.list, function(x) unlist(x$mu1))
all.mu1 <- bind_rows(all.mu1, .id = "ng")
all.mu1.gather <- gather(all.mu1, value = "viables", key = "ng")

all.gather <- list("0.0" = all.mu0.gather, "1.0" = all.mu1.gather)
all.gather <- bind_rows(all.gather, .id = "mu")


ggplot(all.gather, aes(x = mu, y = viables, fill = ng))+
  geom_boxplot(alpha = 0.7)+
  xlab("\U003BC (second gen)")+
  ylab("Viable genomes")+
  scale_y_log10(limits = c(10,100000))+
  scale_x_discrete(labels=c("0.0", "1.0"))+
  scale_fill_manual(values = c("#FFEDA0","#FEB24C","#F03B20"),
                    labels = c("0.0","0.1", "1.0"))+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 19),
        legend.title = element_text(size = 22))+
  labs(fill = "\U003BC (first gen)")
ggsave("figs/next_gen/viables_nextgen.png",
       width = 9, height = 4)

wilcox.test(all.mu0$ng0, all.mu0$ng01, paired = F) #ns
wilcox.test(all.mu0$ng0, all.mu0$ng1, paired = F) #ns
wilcox.test(all.mu1$ng0, all.mu1$ng01, paired = F) #ns
wilcox.test(all.mu1$ng0, all.mu1$ng1, paired = F) #ns

mean(all.mu0$ng0) #518.93
mean(all.mu0$ng01) #533.33

median(all.mu0$ng0) #25
median(all.mu0$ng01) #34


#calculate fitness for each genotype
load("rdata/next_gen_mu0.RData")
load("rdata/next_gen_mu01.RData")
load("rdata/next_gen_mu1.RData")
load("rdata/combo.RData")

calc.fitness.noepi <- function(gen){
  fitness <- as.data.frame(matrix(ncol = 8, nrow = 1))
  colnames(fitness) <- combo
  fitness[1] <- min(c(gen[[1]][1,1],
                      gen[[1]][2,1], 
                      gen[[1]][3,1])) #A1B1C1
  fitness[2] <- min(c(gen[[1]][1,1],
                      gen[[1]][2,1], 
                      gen[[2]][3,1])) #A1B1C2
  fitness[3] <- min(c(gen[[1]][1,1],
                      gen[[2]][2,1], 
                      gen[[1]][3,1])) #A1B2C1
  fitness[4] <- min(c(gen[[1]][1,1],
                      gen[[2]][2,1], 
                      gen[[2]][3,1])) #A1B2C2
  fitness[5] <- min(c(gen[[2]][1,1],
                      gen[[1]][2,1], 
                      gen[[1]][3,1])) #A2B1C1
  fitness[6] <- min(c(gen[[2]][1,1],
                      gen[[1]][2,1], 
                      gen[[2]][3,1])) #A2B1C2
  fitness[7] <- min(c(gen[[2]][1,1],
                      gen[[2]][2,1], 
                      gen[[1]][3,1])) #A2B2C1
  fitness[8] <- min(c(gen[[2]][1,1],
                      gen[[2]][2,1], 
                      gen[[2]][3,1])) #A2B2C2
  return(fitness)
}
mu0.fitness <- map_dfr(next.gen.mu0, calc.fitness.noepi)
mu01.fitness <- map_dfr(next.gen.mu01, calc.fitness.noepi)
mu1.fitness <- map_dfr(next.gen.mu1, calc.fitness.noepi)

ng.fitness.list <- list("ng0" = mu0.fitness, "ng01" = mu01.fitness, "ng1" = mu1.fitness)

save(ng.fitness.list, file = "rdata/nextgen_fitness_list.RData")

#calculate weighted fitness per population
load("rdata/next_gen_counts_list.RData")

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


ng.weighted.list <- map2(ng.fitness.list, ng.counts.list, apply.weighted.fitness)
ng.weighted.list <- map(ng.weighted.list, keep_at, c("mu0","mu1"))
ng.weighted.list <- map(ng.weighted.list, bind_rows, .id = "mu")

ng.weighted.all <- bind_rows(ng.weighted.list, .id = "ng")
ng.weighted.gather <- gather(ng.weighted.all, value = "fit", key = "mu", -ng)

ggplot(ng.weighted.gather, aes(x = mu, y = fit, fill = ng))+
  geom_boxplot(alpha = 0.7)+
  xlab("\U003BC (second gen)")+
  ylab("Weighted fitness")+
  ylim(c(0,1.5))+
  scale_x_discrete(labels=c("0.0", "1.0"))+
  scale_fill_manual(values = c("#FFEDA0","#FEB24C","#F03B20"),
                    labels = c("0.0","0.1", "1.0"))+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 19),
        legend.title = element_text(size = 22))+
  labs(fill = "\U003BC (first gen)")
ggsave("figs/next_gen/fitness_nextgen.png",
       width = 9, height = 4)

wilcox.test(ng.weighted.list$ng0$mu0, ng.weighted.list$ng01$mu0, paired = F) #ns
wilcox.test(ng.weighted.list$ng0$mu0, ng.weighted.list$ng1$mu0, paired = F) #ns
wilcox.test(ng.weighted.list$ng0$mu1, ng.weighted.list$ng01$mu1, paired = F) #ns
wilcox.test(ng.weighted.list$ng0$mu1, ng.weighted.list$ng1$mu1, paired = F) #ns

mean(ng.weighted.list$ng0$mu0) #0.2572398
mean(ng.weighted.list$ng01$mu0) #0.2736658

median(ng.weighted.list$ng0$mu0) #0.1869751
median(ng.weighted.list$ng01$mu0) #0.2098971


#compare parental proportions
load("rdata/parental_matrix_ng01.RData")
load("rdata/parental_matrix_ng1.RData")
parental.matrix.ng01 <- map(parental.matrix.ng01, matrix, nrow  = 1, ncol = 2)
parental.matrix.ng01 <- map(parental.matrix.ng01, as.data.frame)
parental.matrix.ng01 <- bind_rows(parental.matrix.ng01)

ng01.mu0 <- final.gen.ng01$mu0
ng01.mu0 <- ng01.mu0[,c(1,8)]

parental.ng01.1 <- ng01.mu0$A1B1C1[which(parental.matrix.ng01$V1 == TRUE)]
parental.ng01.2 <- ng01.mu0$A2B2C2[which(parental.matrix.ng01$V2 == TRUE)]

novel.ng01.1 <- ng01.mu0$A1B1C1[which(parental.matrix.ng01$V1 == FALSE)]
novel.ng01.2 <- ng01.mu0$A2B2C2[which(parental.matrix.ng01$V2 == FALSE)]

ng01.novel <- sum(novel.ng01.1, novel.ng01.2)/
  (sum(novel.ng01.1, novel.ng01.2,parental.ng01.1, parental.ng01.2))
#0.1391798
ng01.starting.novel <- (length(which(parental.matrix.ng01$V1 == FALSE))+ 
                          length(which(parental.matrix.ng01$V2 == FALSE)))/(length(parental.matrix.ng01$V1) + length(parental.matrix.ng01$V2))
#0.27

parental.matrix.ng1 <- map(parental.matrix.ng1, matrix, nrow  = 1, ncol = 2)
parental.matrix.ng1 <- map(parental.matrix.ng1, as.data.frame)
parental.matrix.ng1 <- bind_rows(parental.matrix.ng1)

ng1.mu0 <- final.gen.ng1$mu0
ng1.mu0 <- ng1.mu0[,c(1,8)]

parental.ng1.1 <- ng1.mu0$A1B1C1[which(parental.matrix.ng1$V1 == TRUE)]
parental.ng1.2 <- ng1.mu0$A2B2C2[which(parental.matrix.ng1$V2 == TRUE)]

novel.ng1.1 <- ng1.mu0$A1B1C1[which(parental.matrix.ng1$V1 == FALSE)]
novel.ng1.2 <- ng1.mu0$A2B2C2[which(parental.matrix.ng1$V2 == FALSE)]

ng1.novel <- sum(novel.ng1.1, novel.ng1.2)/
  (sum(novel.ng1.1, novel.ng1.2,parental.ng1.1, parental.ng1.2))
#0.646214
ng1.starting.novel <- (length(which(parental.matrix.ng1$V1 == FALSE))+ 
                         length(which(parental.matrix.ng1$V2 == FALSE)))/(length(parental.matrix.ng1$V1) + length(parental.matrix.ng1$V2))
#0.73


compare.start.ng01 <- as.data.frame(matrix(ncol = 2, nrow = 2))
colnames(compare.start.ng01) <- c("id", "freq")
compare.start.ng01$id <- c("Parental","Novel")
compare.start.ng01$freq <- c((1- ng01.starting.novel),ng01.starting.novel)

compare.end.ng01 <- as.data.frame(matrix(ncol = 2, nrow = 2))
colnames(compare.end.ng01) <- c("id", "freq")
compare.end.ng01$id <- c("Parental","Novel")
compare.end.ng01$freq <- c((1- ng01.novel),ng01.novel)


ggplot(data = compare.start.ng01, aes(x = "", y = freq, fill = id))+
  geom_col(alpha = 0.9) +
  coord_polar(theta = "y")+
  theme_void()+
  scale_fill_manual(values = c("light blue","#A0B2D8"))+ 
  theme(legend.text = element_text(size = 19),
        legend.title = element_text(size = 22),
        plot.title = element_text(size = 22))
ggsave("figs/next_gen/start_ng01.png",
       width = 4, height = 4, dpi = 300)

ggplot(data = compare.end.ng01, aes(x = "", y = freq, fill = id))+
  geom_col(alpha = 0.9) +
  coord_polar(theta = "y")+
  theme_void()+
  scale_fill_manual(values = c("light blue","#A0B2D8"))+ #"#E8A093"
  theme(legend.text = element_text(size = 19),
        legend.title = element_text(size = 22))
ggsave("figs/next_gen/end_ng01.png",
       width = 4, height = 4, dpi = 300)

compare.start.ng1 <- as.data.frame(matrix(ncol = 2, nrow = 2))
colnames(compare.start.ng1) <- c("id", "freq")
compare.start.ng1$id <- c("Parental","Novel")
compare.start.ng1$freq <- c((1- ng1.starting.novel),ng1.starting.novel)

compare.end.ng1 <- as.data.frame(matrix(ncol = 2, nrow = 2))
colnames(compare.end.ng1) <- c("id", "freq")
compare.end.ng1$id <- c("Parental","Novel")
compare.end.ng1$freq <- c((1- ng1.novel),ng1.novel)

ggplot(data = compare.start.ng1, aes(x = "", y = freq, fill = id))+
  geom_col(alpha = 0.9) +
  coord_polar(theta = "y")+
  theme_void()+
  scale_fill_manual(values = c("light blue","#A0B2D8"))+ #"#E8A093"
  theme(legend.text = element_text(size = 19),
        legend.title = element_text(size = 22))
ggsave("figs/next_gen/start_ng1.png",
       width = 4, height = 4, dpi = 300)

ggplot(data = compare.end.ng1, aes(x = "", y = freq, fill = id))+
  geom_col(alpha = 0.9) +
  coord_polar(theta = "y")+
  theme_void()+
  scale_fill_manual(values = c("light blue","#A0B2D8"))+ #"#E8A093"
  theme(legend.text = element_text(size = 19),
        legend.title = element_text(size = 22))
ggsave("figs/next_gen/end_ng1.png",
       width = 4, height = 4, dpi = 300)

