# finemapping-pipeline

Fine-mapping pipeline using FINEMAP and SuSiE.


# Output descriptions

## SUSIE output

### PHENONAME.SUSIE.cred.bgz
### Contains credible set summaries from SUSIE finemapping for all genome-wide significant regions.

Columns:
- region:  region for which the finemapping was ran. 
- cs: running number for independent credible sets in a region
- cs_log10bf: Log10 bayes factor of comparing the solution of this model (cs  independent credible sets) to cs -1 credible sets.
- cs_avg_r2: Average correlation R2 between variants in the credible set
- cs_min_r2: minimum r2  between variants in the credible set
- cs_size: how many snps does this credible set contain

### PHENONAME.SUSIE.snp.bgz
### Contains variants and their summary statistics for variants in all credible set.

Columns:
- trait: phenotype
- region:  region for which the finemapping was ran.
- v, rsid: variant ids
- chromosome     
-  position        
- allele1 
- allele2 
- maf: minor allele frequency     
- beta    original nbeta 
- se      original se
- p       original p
- mean     posterior mean beta after finemapping 
- sd     posterior standard deviation after finemapping.
- prob    posterior inclusion probability
- cs n credible set within region

## Finemap output

### PHENONAME.FINEMAP.config.bgz
### Summary finemapping variant configurations from FINEMAP method.
Columns:
- trait: phenotype
- region:  region for which the finemapping was ran.
- rank: rank of this configuration within a region
- config: causal variants in this configuration
- prob: probability across all n independent signal configurations
- log10bf: log10 bayes factor for this configuration
- odds: odds of this configuration
- k how many independent signals in this configuration
- prob_norm_k: probability of this configuration within k independent signals solution
- h2 snp heritability of this solution.
- 95% confidence interval limits of snp heritability of this solution.
- mean marginalized shrinkage estimates of the posterior effect size mean
- sd marginalized shrinkage estimates of the posterior effect standard deviation 

### FINEMAP.region.bgz
### Summary statistics on number of independent signals in each region

Columns:
- trait: phenotype
- region:  region for which the finemapping was ran.
- h2g snp heritability of this region
- h2g_sd standard deviation of snp heritability of this region
- h2g_lower95 lower limit of 95% CI for snp heritability
- h2g_upper95 upper limit of 95% CI for snp heritability
- log10bf log bayes factor compared against null (no signals in the region)
- prob_xSNP columns for probabilities of different number of independent signals 
-  expectedvalue: expectation (average) of the number of signals


### PHENOTYPE.FINEMAP.SNP.bgz
### Summary statistics of variants and in to what credible set they may belong to.

Columns:
- trait: phenotype
- region:  region for which the finemapping was ran.
- v: variant
- index: running index
- rsid    variant id
-chromosome      
-position        
-allele1 
-allele2 
-maf     
-beta    
-se      
-z       
-prob    
-log10bf 
- mean marginalized shrinkage estimates of the posterior effect size mean
- sd marginalized shrinkage estimates of the posterior effect standard deviation     

-mean_inclsd_incl:
p       cs      cs2     cs3     cs4     cs5     cs6     cs7     cs8     cs9     cs10

###PHENOTYPE.REGION.credx files
###Posterior inclusion probabilities of SNPs in x number of signals solution

Columns:
-index running index 
-credn variant in credible set n
-probn posterior inclusion probability of variant into credible set n




