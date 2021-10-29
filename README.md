# Statistical fine-mapping pipeline in FinnGen
Fine-mapping pipeline using [FINEMAP](http://finemap.me/) and [SuSiE](https://github.com/stephenslab/susieR).

## Overview
This repository provides a pipeline for statistical fine-mapping in FinnGen. We used two state-of-the-art methods, [FINEMAP](http://finemap.me/) ([Benner, C. et al., 2016](https://academic.oup.com/bioinformatics/article/32/10/1493/1743040); [Benner, C. et al., 2018](https://www.biorxiv.org/content/10.1101/318618v1)) and [SuSiE](https://github.com/stephenslab/susieR) ([Wang, G. et al., 2020](https://www.biorxiv.org/content/10.1101/501114v4)) to fine-map genome-wide significant loci in FinnGen endpoints.

Briefly, there are three main steps:
### 1. Preprocessing
For each genome-wide significant locus (default configuration: P < 5e-8), we define a fine-mapping region by taking a 3 Mb window around a lead variant (and merge regions if they overlap). We preprocess an input GWAS summary statistics into separate files per region for the following steps.

### 2. LD computation
We compute in-sample dosage LD using [LDstore2](http://www.christianbenner.com/) for each fine-mapping region.

### 3. Fine-mapping
With the inputs of summary statistics and in-sample LD from the steps 1-2, we conduct fine-mapping using [FINEMAP](http://finemap.me/) and [SuSiE](https://github.com/stephenslab/susieR) with the maximum number of causal variants in a locus L = 10.

## Usage
Please run [`finemap.wdl`](wdl/finemap.wdl) on a cromwell server with an appropriate json input (e.g., [`finemap_inputs.json`](wdl/finemap_inputs.json)) and a dependent sub workflow ([`finemap_sub.zip`](wdl/finemap_sub.zip)). Please pay attention to an important configuration below (`finemap.ldstore_finemap.susie.min_cs_corr`) controlling reliability of credible sets and variants!

Configurable options include:
- `finemap.sumstats_pattern`: a path to GWAS summary statistics where `{PHENO}` is a magic keyword to be replaced by an actual phenotype name in `finemap.phenolistfile`
- `finemap.phenolistfile`: a list of phenotypes to fine-map
- `finemap.phenotypes`: a path to phenotype file. **[IMPORTANT]** The first column should contains sample IDs that match to bgen sample IDs. If a phenotype is binary, we expect it's coded as **0/1/NA**.
- `finemap.preprocess.rsid_col`: a column name of rsid / variant id. If `finemap.preprocess.set_variant_id` is `true`, this option is disregarded.
- `finemap.preprocess.chromosome_col`: a column name of chromosome
- `finemap.preprocess.position_col`: a column name of position
- `finemap.preprocess.allele1_col`: a column name of reference allele (by default)
- `finemap.preprocess.allele2_col`: a column name of alternative allele (by default)
- `finemap.preprocess.freq_col`: a column name of allele frequency
- `finemap.preprocess.beta_col`: a column name of marginal beta
- `finemap.preprocess.se_col`: a column name of standard error of marginal beta
- `finemap.preprocess.p_col`: a column name of p-value
- `finemap.preprocess.delimiter`: a delimiter of summary statistics. Due to the limitation of non-supported escape characters in json, it supports the following magic keywords:
    - `WHITESPACE`: `\s+`
    - `SINGLE_WHITESPACE`: `'\s'`
    - `TAB`: `'\t'`
    - `SPACE`: `' '`
- `finemap.preprocess.scale_se_by_pval`: an option to scale standard error based on p-value
- `finemap.preprocess.x_chromosome`: an option to include X-chromosome. It assumes females coded as 0/1/2 and males coded as 0/2.
- `finemap.preprocess.set_variant_id`: an option to specify whether to set variant ids as `chr:pos:ref:alt`. If `false`, `finemap.preprocess.rsid_col` is required.
- `finemap.preprocess.max_region_width`: Overlapping regions won't be merged if the resulting region exceeds this width
- `finemap.preprocess.p_threshold`: a p-value threshold to define a fine-mapping region
- `finemap.ldstore_finemap.n_causal_snps`: a maximum number of causal variants per locus
- **`finemap.ldstore_finemap.susie.min_cs_corr`**: **[IMPORTANT]** a minimum pairwise correlation value (`r`) for variants in a credible set for purity filter in SuSiE. In a minority of credible sets, there is a region with a few lead variants in tight LD and low LD to others in a region. In these occasions, if the sum of PIPs of those in tight LD is below 0.95, SuSiE adds low LD variants to get 95% credible sets. However, since those low LD variants are not part of "pure"/reliable credible set, purity filtering filters the whole credible set if minimum r2 between variants is lower than the given threshold. To enable a post-hoc purity filtering, it is set as 0 by default but users are *strongly encouraged* to do a purity filtering based on cs_min_r2 (default original SuSiE filter 0.5) value or low_purity flag. In many occasions, the lead tight LD variants form truly the credible set of variants, and one option, depending on use case, would be post-hoc filtering variants in a credible set by r2 values (which results in < 95% credible set).

- `finemap.ldstore_finemap.combine.good_cred_r2`: good quality credible set minimum r2 between variants. Controls how filtered summary of credible sets and their variants are reported. See output of [filtered credible set report](#PHENONAME.cred.summary.tsv) ,
- `finemap.ldstore_finemap.combine.snp_annot_file`: optional file for arbitrary variant annotations to be included in filtered summary credible set reports. Must be bgzipped and tabixed. Remove the whole row if no annotations are to be used. First four columns must be chr pos ref alt
-`finemap.ldstore_finemap.combine.snp_annot_file_tbi`: tabix index file for the above file. Must be given if annot_file is specified
- `finemap.ldstore_finemap.combine.snp_annot_fields`: comma separated string of column names to include from the annot_file

- `finemap.set_variant_id_map_chr`: comma separated list of chromosome name mappings to perform for summary stat files. Can be useful if bgen file and summary file are coded differently. In output the chr are reported in the original format given in the summary stat file. Example 23=X,24=MT to code 23 to X and 24 to MT
- `finemap.ldstore_finemap.filter_and_summarize.good_cred_r2`: r2 value threshold of minimum pairwise r2 between cs variants to call good CS in summary SUSIE files.
- `finemap.ldstore_finemap.filter_and_summarize.snp_annot_file`: Optional annotation file to annotate variants in SUSIE summarized results. Has to be tbi indexed and have column chr pos ref alt matching those in summary stats and any additional annotation columns
- `finemap.ldstore_finemap.filter_and_summarize.snp_annot_file_tbi`: must be given if annotation file is given.
- `finemap.ldstore_finemap.filter_and_summarize.snp_annotation_fields`: which fields to add from annotation file to variant output
- `finemap.ldstore_finemap.filter_and_summarize.cpu`: number of CPUs to use for filtering and annotation
- `finemap.ldstore_finemap.filter_and_summarize.mem`: amount of RAM to use for filtering and annotation


## Output descriptions

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

#### PHENONAME.SUSIE.cred.summary.tsv
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

#### PHENONAME.SUSIE_99.cred.summary.tsv
The same file as PHENONAME.SUSIE.cred.summary.tsv except 99% credible set computed instead of 95%.


#### PHENONAME.SUSIE_extend.cred.summary.tsv
The same file as PHENONAME.SUSIE.cred.summary.tsv (i.e. 95% CS results) except the corresponding SNPs (see below) contain
variants up to 99% credible set.



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

#### PHENONAME.SUSIE.snp.filter.tsv
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

#### PHENONAME.SUSIE_99.snp.filter.tsv
Same as PHENONAME.SUSIE.snp.filter.tsv but corresponding to 99% CS results.

#### PHENONAME.SUSIE_extend.snp.filter.tsv
Same as PHENONAME.SUSIE.snp.filter.tsv but may contain extra variants on top of 95% CS up to forming 99% CS.
Contains 2 additional columns:

- cs_99 indicator which CS this variant is part of.
- cs_specific_prob_99 PIP in 99% CS solution
Note that in case variant is not part of 95% CS but belongs to 99% CS, cs field=-1 and cs_99 contains the CS id.



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

#### FINEMAP.region.bgz
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

#### PHENOTYPE.REGION.credx files
Posterior inclusion probabilities of SNPs in x number of signals solution

Columns:
- index running index
- credn variant in credible set n
- probn posterior inclusion probability of variant into credible set n

## Authors
Masahiro Kanai (mkanai@broadinstitute.org) with significant inputs from Hilary Finucane, Juha Karjalainen, Mitja Kurki, Sina Rüeger, and Jacob Ulirsch.
