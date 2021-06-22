# Statistical fine-mapping in FinnGen release 6

Two set of results provided for each phenotype with significant results.
SUSIE based pre-filtered results for good quality credible set variants. All subsequent FinnGen analyses (colocalization, autoreporting) are based on these results. Location: gs://finngen-production-library-green/finngen_R6/finngen_R6_analysis_data/finemap/summaries/

Full SUSIE and FINEMAP results in gs://finngen-production-library-green/finngen_R6/finngen_R6_analysis_data/finemap/

See below for full explanation of all provided files.

## Overview
We used two state-of-the-art methods, [FINEMAP](http://finemap.me/) ([Benner, C. et al., 2016](https://academic.oup.com/bioinformatics/article/32/10/1493/1743040); [Benner, C. et al., 2018](https://www.biorxiv.org/content/10.1101/318618v1)) and [SuSiE](https://github.com/stephenslab/susieR) ([Wang, G. et al., 2020](https://www.biorxiv.org/content/10.1101/501114v4)) to fine-map genome-wide significant loci in FinnGen endpoints.

Briefly, there are three main steps:
### 1. Preprocessing
For each genome-wide significant locus (default configuration: P < 5e-8), we define a fine-mapping region by taking a 3 Mb window around a lead variant (and merge regions if they overlap). We preprocess an input GWAS summary statistics into separate files per region for the following steps.

### 2. LD computation
We compute in-sample dosage LD using [LDstore2](http://www.christianbenner.com/) for each fine-mapping region.

### 3. Fine-mapping
With the inputs of summary statistics and in-sample LD from the steps 1-2, we conduct fine-mapping using [FINEMAP](http://finemap.me/) and [SuSiE](https://github.com/stephenslab/susieR) with the maximum number of causal variants in a locus L = 10.

## File descriptions

Two sets of results are provided for each phenotype with significant associations.

### SUSIE filtered summaries
For each phenotype, 2 files are provided. PHENONAME.cred.summary.tsv (see column description below) gives consensed summary of each credible set (independent association) and associated summary statistics. PHENONAME.snp.filter.tsv (see below) file lists all variants forming the good quality credible sets given in PHENONAME.cred.summary.tsv. In case of good_cs=False no SNPs are listed.

Location of files: gs://finngen-production-library-green/finngen_R6/finngen_R6_analysis_data/finemap/summaries/

### Full finemap files of SUSIE and FINEMAP
Location of files: gs://finngen-production-library-green/finngen_R6/finngen_R6_analysis_data/finemap/full/
See below for explanation of SUSIE (PHENONAME.SUSIE.snp.bgz, PHENONAME.SUSIE.cred.bgz) and FINEMAP (PHENONAME.FINEMAP.config.bgz,PHENONAME.FINEMAP.region.bgz,PHENOTYPE.FINEMAP.snp.bgz)


### SuSiE outputs

#### PHENONAME.SUSIE.cred.bgz
Contains credible set summaries from SuSiE fine-mapping for all genome-wide significant regions.

Columns:
- region: region for which the fine-mapping was run.
- cs: running number for independent credible sets in a region
- cs_log10bf: Log10 bayes factor of comparing the solution of this model (cs independent credible sets) to cs -1 credible sets.
- cs_avg_r2: Average correlation R2 between variants in the credible set
- cs_min_r2: minimum r2 between variants in the credible set
- cs_size: how many snps does this credible set contain

#### PHENONAME.cred.summary.tsv
Summary of credible sets where top variant for each CS is included.

Columns:
- trait: phenotype
- region: region for which the fine-mapping was run.
- cs: running number for independent credible sets in a region
- cs_log10bf: Log10 bayes factor of comparing the solution of this model (cs independent credible sets) to cs -1 credible sets
- cs_avg_r2:	Average correlation R2 between variants in the credible set
- cs_min_r2:	minimum r2 between variants in the credible set
- low_purity: boolean (TRUE,FALSE) indicator if the CS is low purity (low min r2)
- cs_size: how many snps does this credible set contain
- good_cs: boolean (TRUE,FALSE) indicator if this CS is considered reliable. IF this is FALSE then top variant reported for the CS will be chosen based on minimum p-value in the credible set, otherwise top variant is chosen by maximum PIP
- cs_id:
- v: top variant (chr:pos:ref:alt)
- p: top variant p-value
- beta: top variant beta
- sd: top variant standard deviation
- prob: overall PIP of the variant in the region
- cs_specific_prob:	PIP of the variant in the current credible set (this and previous are typically almost identical)
- 0..n: configured annotation columns. Typical default most_severe,gene_most_severe giving consequence and gene of top variant


#### PHENONAME.SUSIE.snp.bgz
Contains variant summaries with credible set information.

Columns:
- trait: phenotype
- region: region for which the fine-mapping was run.
- v, rsid: variant ids
- chromosome
- position
- allele1
- allele2
- maf: minor allele frequency
- beta: original marginal beta
- se: original se
- p: original p
- mean: posterior mean beta after fine-mapping
- sd: posterior standard deviation after fine-mapping.
- prob: posterior inclusion probability
- cs: credible set index within region
- lead_r2: r2 value to a lead variant (the one with maximum PIP) in a credible set
- alphax: posterior inclusion probability for the x-th single effect (x := 1..L where L is the number of single effects (causal variants) specified; default: L = 10).

#### PHENONAME.snp.filter.tsv
snps that are part of good quality credible sets as reported in `PHENONAME.cred.summary.tsv` file.

- trait phenotype
- region region for which the fine-mapping was run
- v	variantid (chr:pos:ref:alt)
- cs running credible set id within region
- cs_specific_prob	posterior inclusion probability for this CS.
- chromosome
- position
- allele1
- allele2
- maf
- beta	original association beta
- p	original pvalue
- se	original se
- most_severe most severe consequence of the variant
- gene_most_severe gene corresponding to most severe consequence

### FINEMAP outputs

#### PHENONAME.FINEMAP.config.bgz
Summary fine-mapping variant configurations from FINEMAP method.

Columns:
- trait: phenotype
- region: region for which the fine-mapping was run.
- rank: rank of this configuration within a region
- config: causal variants in this configuration
- prob: probability across all n independent signal configurations
- log10bf: log10 bayes factor for this configuration
- odds: odds of this configuration
- k: how many independent signals in this configuration
- prob_norm_k: probability of this configuration within k independent signals solution
- h2: snp heritability of this solution.
- 95% confidence interval limits of snp heritability of this solution.
- mean: marginalized shrinkage estimates of the posterior effect size mean
- sd: marginalized shrinkage estimates of the posterior effect standard deviation

#### PHENONAME.FINEMAP.region.bgz
Summary statistics on number of independent signals in each region

Columns:
- trait: phenotype
- region: region for which the fine-mapping was run.
- h2g snp: heritability of this region
- h2g_sd: standard deviation of snp heritability of this region
- h2g_lower95: lower limit of 95% CI for snp heritability
- h2g_upper95: upper limit of 95% CI for snp heritability
- log10bf: log bayes factor compared against null (no signals in the region)
- prob_xSNP: columns for probabilities of different number of independent signals
- expectedvalue: expectation (average) of the number of signals

#### PHENOTYPE.FINEMAP.snp.bgz
Summary statistics of variants and into what credible set they may belong to.

Columns:
- trait: phenotype
- region: region for which the fine-mapping was run.
- v: variant
- index: running index
- rsid: variant id
- chromosome
- position
- allele1
- allele2
- maf: minor allele frequency
- beta: original marginal beta
- se: original se
- z: original z
- prob: posterior inclusion probability
- log10bf: log10 bayes factor
- mean: marginalized shrinkage estimates of the posterior effect size mean
- sd: marginalized shrinkage estimates of the posterior effect standard deviation
- mean_incl: conditional estimates of the posterior effect size mean
- sd_incl: conditional estimates of the posterior effect size standard deviation
- p: original p
- csx: credible set index for given number of causal variants x
