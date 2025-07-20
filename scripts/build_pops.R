library(tidyverse)
library(pegas)
library(here)

#generate fitness distribution
dms.dist <- read_csv("dms_dist.csv")
colnames(dms.dist) <- "x"
dms.dist <- na.omit(dms.dist)
dms.dist[which(dms.dist$x <= 0),] <- 0.00001

ggplot(data = dms.dist, aes(x = x))+
  geom_histogram(bins = 30, fill = "grey", color = "white")+
  scale_x_continuous(limits = c(0, 1.7), oob = scales::oob_keep)+
  ylab("Count")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        plot.title = element_text(size = 22))
ggsave("runs/run_051525/figs/fitness_dist.png", width = 5, height = 4, dpi = 300)

#generate fitness distributions under immune pressure
viable.dist <- dms.dist[which(dms.dist$x > 0.00001),]
set.seed(1)
dist.5 <- as.data.frame(c(rep(0.00001,950),sample(viable.dist$x,50))) #5% resistance
colnames(dist.5) <- "x"
set.seed(1)
dist.15 <- as.data.frame(c(rep(0.00001,850),sample(viable.dist$x,150))) #15% resistance
colnames(dist.15) <- "x"
set.seed(1)
dist.25 <- as.data.frame(c(rep(0.00001,750),sample(viable.dist$x,250))) #25% resistance
colnames(dist.25) <- "x"
set.seed(1)
dist.35 <- as.data.frame(c(rep(0.00001,650),sample(viable.dist$x,350))) #35% resistance
colnames(dist.35) <- "x"

#15% resistance
ggplot(data = dist.15, aes(x = x))+
  geom_histogram(bins = 30, fill = "grey", color = "white")+
  scale_x_continuous(limits = c(0, 1.7), oob = scales::oob_keep)+
  ylab("Count")+
  scale_y_log10(limits = c(1,1000))+
  ggtitle("15% resistance")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        plot.title = element_text(size = 22))
ggsave("runs/run_051525/figs/15percent_dist.png", width = 5, height = 4, dpi = 300)

#generate relatedness spectra
pb2 <- read.FASTA("clustalo_PB2_alignment.fa")
h1n1.pb2.seq <- pb2[grep("H1N1", names(pb2))]
h3n2.pb2.seq <- pb2[grep("H3N2", names(pb2))]

seq.compare <- function(h3,h1){
  concat <- list()
  div <- list()
  for(i in seq_along(h3)){
    concat[[i]] <- c(h3[i],h1)
    div[[i]] <- nuc.div(concat[[i]])
  }
  div <- unlist(div)
  return(div)
}

h1n1.h3n2.pb2.div <- list()
for(i in seq_along(h1n1.pb2.seq)){
  h1n1.h3n2.pb2.div[[i]] <- seq.compare(h3n2.pb2.seq,h1n1.pb2.seq[i])
}
mean(unlist(h1n1.h3n2.pb2.div))#mean = 0.1690215 
sd(unlist(h1n1.h3n2.pb2.div)) #sd = 0.001427752

h1n1.pb2.div <- nuc.div(h1n1.pb2.seq, variance = T) #mean = 0.007161852 #sd = 0.003548207
h3n2.pb2.div <- nuc.div(h3n2.pb2.seq, variance = T) #mean = 0.004649173 #sd = 0.002346407

set.seed(1)
h1.diff <- rnorm(100, mean = h1n1.pb2.div[1], sd = sqrt(h1n1.pb2.div[2]))
set.seed(1)
h3.diff <- rnorm(100, mean = h3n2.pb2.div[1], sd = sqrt(h3n2.pb2.div[2]))
set.seed(1)
h1h3.diff <- rnorm(200, mean = mean(unlist(h1n1.h3n2.pb2.div)), sd = sd(unlist(h1n1.h3n2.pb2.div)))

diffs <- c(h1.diff, h3.diff, h1h3.diff)
diffs.norm <- as.data.frame((diffs - min(diffs))/(max(diffs) - min(diffs)))
colnames(diffs.norm) <- "r"

ggplot(data = diffs.norm, aes(x = r))+
  geom_histogram(bins = 50, fill = "grey", color = "white")+
  scale_x_continuous(limits = c(0, 1), oob = scales::oob_keep)+
  xlab("∆r")+
  ylab("Count")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        plot.title = element_text(size = 22))
ggsave("runs/run_051525/figs/r_hist.png",width = 5, height = 4, dpi = 300)

#build genotypes
build.abc <- function(fit, diff){
  gen.mat <- list("i1" = NA, "i2" = NA)
  gen.mat <- map(gen.mat, function(x) x <- as.data.frame(matrix(data = NA, nrow = 3, ncol = 2,
                                                                dimnames = list(c("A","B","C"),
                                                                                c("x","r")))))
  gen.mat$i1$x <- sample(fit$x, 3, replace = F)
  gen.mat$i2$x <- sample(fit$x, 3, replace = F)
  gen.mat$i1$r <- 0
  gen.mat$i2$r <- sample(diff$r, 1, replace = F)
  return(gen.mat)
}
build.ab <- function(fit, fit.imm,diff){
  gen.mat <- list("i1" = NA, "i2" = NA)
  gen.mat <- map(gen.mat, function(x) x <- as.data.frame(matrix(data = NA, nrow = 3, ncol = 2,
                                                                dimnames = list(c("A","B","C"),
                                                                                c("x","r")))))
  gen.mat$i1$x[c(1:2)] <- sample(fit.imm$x, 2, replace = F)
  gen.mat$i2$x[c(1:2)] <- sample(fit.imm$x, 2, replace = F)
  gen.mat$i1$x[3] <- sample(fit$x, 1, replace = F)
  gen.mat$i2$x[3] <- sample(fit$x, 1, replace = F)
  gen.mat$i1$r <- 0
  gen.mat$i2$r <- sample(diff$r, 1, replace = F)
  return(gen.mat)
}
build.c <- function(fit, fit.imm,diff){
  gen.mat <- list("i1" = NA, "i2" = NA)
  gen.mat <- map(gen.mat, function(x) x <- as.data.frame(matrix(data = NA, nrow = 3, ncol = 2,
                                                                dimnames = list(c("A","B","C"),
                                                                                c("x","r")))))
  gen.mat$i1$x[c(1:2)] <- sample(fit$x, 2, replace = F)
  gen.mat$i2$x[c(1:2)] <- sample(fit$x, 2, replace = F)
  gen.mat$i1$x[3] <- sample(fit.imm$x, 1, replace = F)
  gen.mat$i2$x[3] <- sample(fit.imm$x, 1, replace = F)
  gen.mat$i1$r <- 0
  gen.mat$i2$r <- sample(diff$r, 1, replace = F)
  return(gen.mat)
}

set.seed(1)
gen.lists <- replicate(100, build.abc(dms.dist, diffs.norm), simplify = F)
save(gen.lists, file = "runs/run_051525/rdata/gen_lists.RData")

#build unbalanced genotypes
build.abc.unb <- function(fit, diff){
  gen.mat <- list("i1" = NA, "i2" = NA)
  gen.mat <- map(gen.mat, function(x) x <- as.data.frame(matrix(data = NA, nrow = 3, ncol = 2,
                                                                dimnames = list(c("A","B","C"),
                                                                                c("x","r")))))
  gen.mat$i1$x <- sample(fit$x, 3, replace = F)
  gen.mat$i2$x <- sample(fit$x, 3, replace = F)
  gen.mat$i1$r[1] <- 0
  gen.mat$i1$r[c(2:3)] <- sample(diff$r, 2, replace = F)
  gen.mat$i2$r <- sample(diff$r, 3, replace = F)
  return(gen.mat)
}

set.seed(1)
gen.lists.unbalanced <- replicate(100, build.abc.unb(dms.dist, diffs.norm), simplify = F)
save(gen.lists.unbalanced, file = "runs/run_051525/rdata/gen_lists_unbalanced.RData")

#different resistance levels
#5% resistance
set.seed(2)
gen.list.ab.5 <- replicate(100, build.ab(dms.dist,dist.5,diffs.norm), simplify = F)
set.seed(3)
gen.list.c.5 <- replicate(100, build.c(dms.dist,dist.5,diffs.norm), simplify = F)
gen.lists.5 <- list("abc" = gen.lists, "ab" = gen.list.ab.5,"c" = gen.list.c.5)
save(gen.lists.5, file = "runs/run_051525/rdata/gen_lists_5.RData")

#15% resistance
set.seed(2)
gen.list.ab.15 <- replicate(100, build.ab(dms.dist,dist.15,diffs.norm), simplify = F)
set.seed(3)
gen.list.c.15 <- replicate(100, build.c(dms.dist,dist.15,diffs.norm), simplify = F)
gen.lists.15 <- list("abc" = gen.lists, "ab" = gen.list.ab.15,"c" = gen.list.c.15)
save(gen.lists.15, file = "runs/run_051525/rdata/gen_lists_15.RData")

#25% resistance
set.seed(2)
gen.list.ab.25 <- replicate(100, build.ab(dms.dist,dist.25,diffs.norm), simplify = F)
set.seed(3)
gen.list.c.25 <- replicate(100, build.c(dms.dist,dist.25,diffs.norm), simplify = F)
gen.lists.25 <- list("abc" = gen.lists, "ab" = gen.list.ab.25,"c" = gen.list.c.25)
save(gen.lists.25, file = "runs/run_051525/rdata/gen_lists_25.RData")

#35% resistance
set.seed(2)
gen.list.ab.35 <- replicate(100, build.ab(dms.dist,dist.35,diffs.norm), simplify = F)
set.seed(3)
gen.list.c.35 <- replicate(100, build.c(dms.dist,dist.35,diffs.norm), simplify = F)
gen.lists.35 <- list("abc" = gen.lists, "ab" = gen.list.ab.35,"c" = gen.list.c.35)
save(gen.lists.35, file = "runs/run_051525/rdata/gen_lists_35.RData")

#generate list of possible genotype combinations
combo <- as.data.frame(matrix(data = c("A1","A2","B1","B2","C1","C2"), ncol = 3, nrow = 2))
colnames(combo) <- c("A","B","C")
combo <- expand(combo,A,B,C)
combo <- split(combo,1:nrow(combo))
combo <- unlist(lapply(combo, paste0, collapse = ""))
save(combo,file = "runs/run_051525/rdata/combo.RData")

