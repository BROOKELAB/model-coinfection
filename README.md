# model-coinfection
Code for simulations published in "The fitness consequences of coinfection and reassortment among segmented viruses depend upon viral genetic structure" (Farjo & Brooke, 2025)

**Scripts:**

Generating coinfection pairs:
  - *build_pops.R*: build viruses from DMS distribution and sequence diversity data
  - *build_next_gen.R*: select coinfection partners from a post-replication pool for a second generation of infection
  
Single-cycle replication and packaging runs:
  - *run_noepistasis.R*: in the absence of epistasis
  - *run_epistasis.R*: with epistasis
  - *run_positive.R*: with positive dominance dynamics
  - *run_negative.R*: with negative dominance dynamics
  - *run_averaged.R*: with replication propensities determined by averaged per-segment fitness

Data analysis:
  - *viable_genomes.R*: calculate total complete genomes per replication cycle
  - *genome_fitness.R*: calculate weighted population fitness
  - *max_fitness.R*: calculate maximum genotypic fitness per population
  - *resistance_levels.R*: calculate total genomes and population fitness under immune pressure
  - *generation_change.R*: calculate changes in allele frequencies post-replication
  - *compare_next_gen.R*: calculate total genomes and population fitness after a second generation of infection

**Data:**

*clustalo_PB2_alignment.fa*: FASTA file containing PB2 sequences of co-circulating influenza viruses during the 2024-2025 flu season

*dms_dist.csv*: influenza fitness data derived from deep mutational scanning (adapted from [Liu et al (2024) *Virus Evolution*](https://academic.oup.com/ve/article/10/1/veae046/7696177))

*gen_lists.RData*: per-segment fitness and phylogenetic distance for coinfection partners
  - *gen_lists_5.RData*: at 5% immune resistance
  - *gen_lists_15.RData*: at 15% immune resistance
  - *gen_lists_25.RData*: at 25% immune resistance
  - *gen_lists_35.RData*: at 35% immune resistance
  - *resistance_gen_list.RData*: all immune resistance lists
  - *next_gen_mu0.RData*: next generation, pulled from μ = 0.0 pool
  - *next_gen_mu0_epi.RData*: next generation, pulled from μ = 0.0 pool + epistasis
  - *next_gen_mu01.RData*: next generation, pulled from μ = 0.1 pool
  - *next_gen_mu01_epi.RData*: next generation, pulled from μ = 0.1 pool + epistasis
  - *next_gen_mu1.RData*: next generation, pulled from μ = 1.0 pool
  - *next_gen_mu1_epi.RData*: next generation, pulled from μ = 1.0 pool + epistasis

*combo.RData*: every possible parental and reassortant genotype

*fitness_list.RData*: per-genotype fitness calculations
  - *resistance_fitness_list.RData*: fitness under different levels of immune resistance
  - *next_gen_fitness_list.RData*: fitness of genotypes selected for a second round of infection

*viables_list.RData*: total viable genomes for each condition (no epistasis, epistasis, positive dominance, negative dominance, averaged fitness)
  - *resistance_viables_list.RData*: total viable genomes under different levels of population immune resistance
  - *next_gen_viables_list.RData*: total viable genomes for a second round of infection

*counts_list.RData*: final counts of each genotype post-replication
  - *resistance_counts_list.RData*: final counts under immune pressure
  - *next_gen_counts_list.RData*: final counts for a second round of infection

*parental_matrix_ng01.RData*: proportion of parental vs novel genotypes selected from the μ = 0.1 pool for the second generation of replication 

*parental_matrix_ng1.RData*: proportion of parental vs novel genotypes selected from the μ = 1.0 pool for the second generation of replication 



  


