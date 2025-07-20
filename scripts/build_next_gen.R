library(tidyverse)
library(here)

load("rdata/counts_list.RData")
load("rdata/gen_lists.RData")

mu0.gen <- counts.list$final.gen.noepi$mu0
mu01.gen <- counts.list$final.gen.noepi$mu01
mu1.gen <- counts.list$final.gen.noepi$mu1

mu0.gen.epi <- counts.list$final.gen.epi$mu0
mu01.gen.epi <- counts.list$final.gen.epi$mu01
mu1.gen.epi <- counts.list$final.gen.epi$mu1

split.genes <- function(select, gen){
  split.s <- list()
  for(i in seq_along(select)){
    split.s[[i]] <- sapply(seq(from=1, to=nchar(select[[i]]), by=2), function(x) substr(select[[i]], x, x+1))
  }
  new.genes <- list("i1" = NA, "i2" = NA)
  new.genes$i1 <- matrix(data = NA, nrow = 3, ncol = 2)
  rownames(new.genes$i1) <- c("A","B","C")
  colnames(new.genes$i1) <- c("x","r")
  new.genes$i2 <- matrix(data = NA, nrow = 3, ncol = 2)
  rownames(new.genes$i2) <- c("A","B","C")
  colnames(new.genes$i2) <- c("x","r")
  new.genes <- map(new.genes, as.data.frame)
  
  for(i in seq_along(rownames(new.genes$i1))){
    if(!is.na(str_extract(split.s[[1]][i],"1"))){
      new.genes$i1[i,] <- gen$i1[i,]
    } else {
      new.genes$i1[i,] <- gen$i2[i,]
    }
    if(!is.na(str_extract(split.s[[2]][i],"1"))){
      new.genes$i2[i,] <- gen$i1[i,]
    } else {
      new.genes$i2[i,] <- gen$i2[i,]
    }
  }
  return(new.genes)
}

select.gens <- function(counts, gen){
  gen.reps <- list()
  gen.select <- list()
  for(i in seq_along(rownames(counts))){
    gen.reps[[i]] <- rep(names(as.list(counts[i,])), as.list(counts[i,]))
    gen.select[[i]] <- sample(gen.reps[[i]], 2, replace = F)
  }
  next.gen <- map2(gen.select, gen, split.genes)
  for(i in seq_along(next.gen)){
    names(next.gen[[i]]) <- gen.select[[i]]
  }
  return(next.gen)
}

set.seed(1)
next.gen.mu0 <- select.gens(mu0.gen, gen.lists)
save(next.gen.mu0, file = "rdata/next_gen_mu0.RData")

set.seed(1)
next.gen.mu01 <- select.gens(mu01.gen, gen.lists)
save(next.gen.mu01, file = "rdata/next_gen_mu01.RData")

set.seed(1)
next.gen.mu1 <- select.gens(mu1.gen, gen.lists)
save(next.gen.mu1, file = "rdata/next_gen_mu1.RData")

set.seed(1)
next.gen.mu0.epi <- select.gens(mu0.gen.epi, gen.lists)
save(next.gen.mu0.epi, file = "rdata/next_gen_mu0_epi.RData")

set.seed(1)
next.gen.mu01.epi <- select.gens(mu01.gen.epi, gen.lists)
save(next.gen.mu01.epi, file = "rdata/next_gen_mu01_epi.RData")

set.seed(1)
next.gen.mu1.epi <- select.gens(mu1.gen.epi, gen.lists)
save(next.gen.mu1.epi, file = "rdata/next_gen_mu1_epi.RData")

#classify based on parental/novel genotype
names.ng01 <- map(next.gen.mu01, names)
parental.matrix.ng01 <- list()
for(i in seq_along(names.ng01)){
  parental.matrix.ng01[[i]] <- grepl("A1B1C1|A2B2C2",names.ng01[[i]])
}
save(parental.matrix.ng01, file = "rdata/parental_matrix_ng01.RData")

names.ng1 <- map(next.gen.mu1, names)
parental.matrix.ng1 <- list()
for(i in seq_along(names.ng1)){
  parental.matrix.ng1[[i]] <- grepl("A1B1C1|A2B2C2",names.ng1[[i]])
}
save(parental.matrix.ng1, file = "rdata/parental_matrix_ng1.RData")

