library(GillespieSSA)
library(purrr)
library(here)

#load in starting genotype data
load("rdata/gen_lists.RData")

#define general parameters
pop.init <- c(5,5) #number of initiating genomes from genotypes 1 and 2
mu.list <- c(0,0.1,0.5,0.9,1) #possible mixing rates
tf <- 8 #final time

#define replication parameters
build.rep.parms <- function(gen.list){
  rep.parms <- list()
  for(i in seq_along(gen.list)){
    rep.parms[[i]] <- c(x.a1 = gen.list[[i]]$i1[1,1],
                        x.a2 = gen.list[[i]]$i2[1,1],
                        x.b1 = gen.list[[i]]$i1[2,1],
                        x.b2 = gen.list[[i]]$i2[2,1],
                        x.c1 = gen.list[[i]]$i1[3,1],
                        x.c2 = gen.list[[i]]$i2[3,1],
                        r.a1 = gen.list[[i]]$i1[1,2],
                        r.a2 = gen.list[[i]]$i2[1,2],
                        r.b1 = gen.list[[i]]$i1[2,2],
                        r.b2 = gen.list[[i]]$i2[2,2],
                        r.c1 = gen.list[[i]]$i1[3,2],
                        r.c2 = gen.list[[i]]$i2[3,2])
  }
  return(rep.parms)
}

rep.parms <- build.rep.parms(gen.lists)

#initial segment counts
segment.counts <- rep(c(pop.init[1], pop.init[2]), 3)
names(segment.counts) <- c("A1","A2","B1","B2","C1","C2")

segments.list <- list("mu0" = segment.counts, "mu01" = segment.counts,
                      "mu05" = segment.counts, "mu09" = segment.counts,
                      "mu1" = segment.counts) 

segment.calc.1 <- function(segments, mu){
  segments[c(1,3,5)] <- segments[c(1,3,5)]/(1+mu)
  segments[c(2,4,6)] <- segments[c(2,4,6)] - segments[c(1,3,5)]
  return(segments)
}
segment.calc.2 <- function(segments, mu){
  segments[c(2,4,6)] <- segments[c(2,4,6)]/(1+mu)
  segments[c(1,3,5)] <- segments[c(1,3,5)] - segments[c(2,4,6)]
  return(segments)
}

pop1.segments <- map2(segments.list, mu.list, segment.calc.1)
pop2.segments <- map2(segments.list, mu.list, segment.calc.2)

#replication state change matrix
rep.effects <- as.data.frame(matrix(nrow = 6, ncol = 6))
colnames(rep.effects) <- c("A1","A2","B1","B2","C1","C2") 
rownames(rep.effects) <- c("plus.a1","plus.a2","plus.b1","plus.b2","plus.c1","plus.c2")
rep.effects[1,] <- c(1,0,0,0,0,0)
rep.effects[2,] <- c(0,1,0,0,0,0)
rep.effects[3,] <- c(0,0,1,0,0,0)
rep.effects[4,] <- c(0,0,0,1,0,0)
rep.effects[5,] <- c(0,0,0,0,1,0)
rep.effects[6,] <- c(0,0,0,0,0,1)
rep.effects <- t(rep.effects)

#replication propensities
wA1 <- "(x.a1)"
wA2 <- "(x.a2)"
wB1 <- "(x.b1)"
wB2 <- "(x.b2)"
wC1 <- "(x.c1)"
wC2 <- "(x.c2)"

num <- paste("(",wA1,"* A1 +",wA2,"* A2 +",wB1,"* B1 +", wB2,"* B2 +", wC1,"* C1 +",wC2,"* C2)")
denom <- "(A1 + A2 + B1 + B2 + C1 + C2)"

fitness <- paste("(",num,"/",denom,")")

rep.props <- c(
  paste(fitness,"* A1"), #A1
  paste(fitness,"* A2"), #A2
  paste(fitness,"* B1"), #B1
  paste(fitness,"* B2"), #B2
  paste(fitness,"* C1"), #C1
  paste(fitness,"* C2") #C2
)

#simulations
replication.simulation <- function(parms){
  out <- ssa(
    x0 = segments,
    a = rep.props,
    nu = rep.effects,
    parms = parms,
    tf = tf,
    method = ssa.otl(epsilon = 0.01),
    verbose = FALSE
  ) 
  return(out)
}

rep.scheme <- function(parms.scheme){
  rep.data <- map(parms.scheme, replication.simulation)
  rep.data <- map(rep.data, function(x) as.data.frame(x$data))
  return(rep.data)
}

set.seed(1)
pop1.rep <- list("mu0" = NA, "mu01" = NA, "mu05" = NA, "mu09" = NA, "mu1" = NA)
for(i in seq_along(pop1.segments)){
  segments <- pop1.segments[[i]]
  pop1.rep[[i]] <- rep.scheme(rep.parms)
}

set.seed(1)
pop2.rep <- list("mu0" = NA, "mu01" = NA, "mu05" = NA, "mu09" = NA, "mu1" = NA)
for(i in seq_along(pop2.segments)){
  segments <- pop2.segments[[i]]
  pop2.rep[[i]] <- rep.scheme(rep.parms)
}

rep.stats.noepi.avg <- list("pop1" = pop1.rep,"pop2" = pop2.rep)
save(rep.stats.noepi.avg, file = "rdata/rep_stats_noepi_avg.RData")

#calculate total viable genomes
sum.mins <- function(rep.run){
  viable.list <- list()
  for(i in seq_along(rownames(rep.run))){
    viable.list[[i]] <- min((rep.run[i,]$A1 + rep.run[i,]$A2),
                            (rep.run[i,]$B1 + rep.run[i,]$B2),
                            (rep.run[i,]$C1 + rep.run[i,]$C2))
  }
  return(unlist(viable.list))
}

calc.viables <- function(rep.out){
  for(i in seq_along(rep.out)){
    rep.out[[i]] <- map_dfr(rep.out[[i]],tail,1)
  }
  viables <- map(rep.out, sum.mins)
  return(viables)
}

pop1.viables <- calc.viables(pop1.rep)
pop2.viables <- calc.viables(pop2.rep)
subpop.viables.noepi.avg <- list("pop1" = pop1.viables, "pop2" = pop2.viables)
save(subpop.viables.noepi.avg, file = "rdata/subpop_viables_noepi_avg.RData")

viables.noepi.avg <- list("mu0"=NA,"mu01"=NA,"mu05"=NA,"mu09"=NA,"mu1"=NA)
viables.noepi.avg <- map2(pop1.viables, pop2.viables, function(x,y) x+y)

save(viables.noepi.avg, file = "rdata/viables_noepi_avg.RData")

#genome packaging based off the final cycle of replication
load("rdata/rep_stats_noepi_avg.RData")

#define general parameters
tf <- 1000

#define initial segment and genotype counts
genotype.counts <- rep(0,8)
names(genotype.counts) <- c("A1B1C1","A1B1C2","A1B2C1","A1B2C2",
                            "A2B1C1","A2B1C2","A2B2C1","A2B2C2")

build.counts <- function(repx){
  repx.counts <- list()
  for(i in seq_along(repx)){
    repx.counts[[i]] <- c(A1 = repx[[i]][dim(repx[[i]])[1],]$A1,
                          A2 = repx[[i]][dim(repx[[i]])[1],]$A2,
                          B1 = repx[[i]][dim(repx[[i]])[1],]$B1,
                          B2 = repx[[i]][dim(repx[[i]])[1],]$B2,
                          C1 = repx[[i]][dim(repx[[i]])[1],]$C1,
                          C2 = repx[[i]][dim(repx[[i]])[1],]$C2,
                          genotype.counts)
  }
  return(repx.counts)
}

build.mu.counts <- function(counts){
  mu.counts <- list("mu0" = NA, "mu01" = NA, "mu05" = NA, "mu09" = NA, "mu1" = NA)
  for(i in seq_along(counts)){
    mu.counts[[i]] <- build.counts(counts[[i]])
  }
  return(mu.counts)
}

pop1.counts <- build.mu.counts(rep.stats.noepi.avg$pop1)
pop2.counts <- build.mu.counts(rep.stats.noepi.avg$pop2)

#packaging state change matrix
pack.effects <- as.data.frame(matrix(nrow = 8, ncol = 14))
colnames(pack.effects) <- c("A1","A2","B1","B2","C1","C2",
                            "A1B1C1","A1B1C2","A1B2C1","A1B2C2",
                            "A2B1C1","A2B1C2","A2B2C1","A2B2C2") 
rownames(pack.effects) <- c("plus.a1b1c1","plus.a1b1c2","plus.a1b2c1","plus.a1b2c2",
                            "plus.a2b1c1","plus.a2b1c2","plus.a2b2c1","plus.a2b2c2")
pack.effects[1,] <- c(-1,0,-1,0,-1,0,1,0,0,0,0,0,0,0)
pack.effects[2,] <- c(-1,0,-1,0,0,-1,0,1,0,0,0,0,0,0)
pack.effects[3,] <- c(-1,0,0,-1,-1,0,0,0,1,0,0,0,0,0)
pack.effects[4,] <- c(-1,0,0,-1,0,-1,0,0,0,1,0,0,0,0)
pack.effects[5,] <- c(0,-1,-1,0,-1,0,0,0,0,0,1,0,0,0)
pack.effects[6,] <- c(0,-1,-1,0,0,-1,0,0,0,0,0,1,0,0)
pack.effects[7,] <- c(0,-1,0,-1,-1,0,0,0,0,0,0,0,1,0)
pack.effects[8,] <- c(0,-1,0,-1,0,-1,0,0,0,0,0,0,0,1)

pack.effects <- t(pack.effects)

#packaging propensities
pack.props <- c(
  "(A1/(A1+A2)) * (B1/(B1+B2)) * (C1/(C1+C2))", #A1B1C1
  "(A1/(A1+A2)) * (B1/(B1+B2)) * (C2/(C1+C2))", #A1B1C2
  "(A1/(A1+A2)) * (B2/(B1+B2)) * (C1/(C1+C2))", #A1B2C1
  "(A1/(A1+A2)) * (B2/(B1+B2)) * (C2/(C1+C2))", #A1B2C2
  "(A2/(A1+A2)) * (B1/(B1+B2)) * (C1/(C1+C2))", #A2B1C1
  "(A2/(A1+A2)) * (B1/(B1+B2)) * (C2/(C1+C2))", #A2B1C2
  "(A2/(A1+A2)) * (B2/(B1+B2)) * (C1/(C1+C2))", #A2B2C1
  "(A2/(A1+A2)) * (B2/(B1+B2)) * (C2/(C1+C2))"  #A2B2C2
)

#simulations
packaging.simulation <- function(counts){
  out <- list()
  for(i in seq_along(counts)){
    out[[i]] <- ssa(
      x0 = counts[[i]],
      a = pack.props,
      nu = pack.effects,
      tf = tf,
      method = ssa.d(),
      verbose = FALSE
    )
  }
  out <- map(out, function(x) as.data.frame(x$data))
  return(out)
}

set.seed(1)
pop1.pack <- list("mu0" = NA, "mu01" = NA, "mu05" = NA, "mu09" = NA, "mu1" = NA)
for(i in seq_along(pop1.counts)){
  counts <- pop1.counts[[i]]
  pop1.pack[[i]] <- packaging.simulation(counts)
}

set.seed(1)
pop2.pack <- list("mu0" = NA, "mu01" = NA, "mu05" = NA, "mu09" = NA, "mu1" = NA)
for(i in seq_along(pop2.counts)){
  counts <- pop2.counts[[i]]
  pop2.pack[[i]] <- packaging.simulation(counts)
}

pack.stats.noepi.avg <- list("pop1" = pop1.pack, "pop2" = pop2.pack)
save(pack.stats.noepi.avg, file = "rdata/pack_stats_noepi_avg.RData")

#final genotype counts
tab.genotypes <- function(stats.mu){
  finals <- map_dfr(stats.mu, tail,1)
  finals <- finals[,-c(1:7)]
  return(finals)
}

pop1.finals <- map(pack.stats.noepi.avg$pop1, tab.genotypes)
pop2.finals <- map(pack.stats.noepi.avg$pop2, tab.genotypes) 

final.gen.noepi.avg <- map2(pop1.finals, pop2.finals, function(x,y) x+y)

save(final.gen.noepi.avg, file = "rdata/final_gen_noepi_avg.RData")

