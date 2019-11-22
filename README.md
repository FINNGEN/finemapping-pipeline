# Autoreporting tool output readme

This readme contains information about the autoreporting tool results for FinnGen summary statistics.

<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Autoreporting tool output readme](#autoreporting-tool-output-readme)
	- [Contents](#contents)
	- [Used Parameters](#used-parameters)
	- [Detailed File Descriptions](#detailed-file-descriptions)
		- [Group Reports (.top.out)](#group-reports-topout)
		- [Variant Reports (report.out)](#variant-reports-reportout)
		- [Annotated Variants (.annotate.out)](#annotated-variants-annotateout)
		- [Filtered Variant Reports (.fetch.out)](#filtered-variant-reports-fetchout)
- [variant| variant identifier.](#variant-variant-identifier)

<!-- /TOC -->


## Contents

Four files are given for each phenotype:
- Group Reports (.top.out)
- Full Variant Reports (.report.out)
- Annotated Variant Reports (.annotate.out)
- Filtered Variant Reports (.fetch.out)

This folder contains the outputs from running the automatic reporting tool for all phenotypes.

The files are further divided into 4 folders.
The group results are in the folder ```group_reports```. These reports contain an overview of what GWAS Catalog phenotypes if any the variants in the groups matched with, the lead variant p-value, enrichment factor, most severe consequence etc.
The full reports are in the folder ```full_reports```. These contain the complete information for all of the genome-wide significant variants. Information includes annotation information from gnoMAD, finngen annotations, GWAS Catalog matched phenotypes, matched variant information, and group-related information.
The annotated variant reports are in the folder  ```annotated```. These contain only the gws variants and their annotation information from gnoMAD and finngen annotations. The full reports contain all of the information that is available here.
The filtered variant reports are in the folder ```filtered```. These reports contain only the variant information, as well as the group information. The full reports contain all of this information.

All of the reports are tab-separated value (tsv) files, and can be opened using e.g. Libreoffice Calc, or your favourite data analysis environment (R, Python etc).

A recommended way to use these files is:
For a phenotype 'PHENOTYPE',
1) Look if the phenotype has any interesting signals from ```group_reports/PHENOTYPE.gz.top.out```.
2) In case you want a closer look into some of the variants in that phenotype, look up those variants from ```full_reports/PHENOTYPE.gz.report.out```.

The files in folders ```annotated``` and ```filtered``` are of no real interest to anyone other than the tool developers.

## Used Parameters
The parameters that were used for these results are described here.

Parameter name | Value | Description
--- | --- | ---
significance threshold | 5e-8 | significance threshold for variants.
grouping method | ld | Variant grouping method. Plink's ld-clumping for 'ld' and simple location-based clumping for 'simple'.
alternate significance threshold | 0.01 | significance threshold for grouping. Variants in LD with the lead variant of a group can be included into a group if their p-value is smaller than this, and they are closer to the lead variant than the range threshold.
grouping locus width | 1.5 Mb | Variants that are at most this far from a lead variant can be grouped with it, provided sufficiently low p-value and high LD.
Ignore genomic region | 6:23000000-38000000 | Ignore this genomic region from the analysis. In these results, the HLA region was ignored.
grouping ld threshold | 0.05 | r^2 threshold for ld clumping.
GWAS Catalog p-value threshold | e5-8 | GWAS Catalog associations were downloaded if they had a p-value lower than this threshold.

Annotations were acquired from gnoMAD. Finngen Release 4 annotations were used.

## Detailed File Descriptions

### Group Reports (.top.out)

These files contain information about the groups. The used columns are detailed in the following table:
Column name  |  Description
--- | ---
locus_id | The group identifier, containing the chromosome, position, reference and alternate alleles of the group lead variant. This is useful if you want to e.g. filter only this group in the variant report.
chr | group chromosome. All of the variants in this group reside in this chromosome.
start | The smallest position coordinate in the group.
end | The largest position coordinate in the group
enrichment | The Finnish Enrichment value of the group lead variant.
most_severe_gene | The gene in which the most severe consequence for the group lead variant is in.
most_severe_consequence | The most severe consequence for the lead variant.
lead_pval | p-value of the group lead variant.
matching_pheno_gwas_catalog_hits | Currently empty. In case interesting phenotypes were supplied to the autoreporting tool, this column would contain them in case one or multiple variants in the group had a match from GWAS Catalog.
other_gwas_hits | This column contains all of the phenotypes, which have been associated with one or more variants in the group, and reported in GWAS Catalog. If this column is empty, it is possible that the variant is novel.

### Variant Reports (report.out)

These files contain the full reports of the gws variants. They contain the variant basic information, such as p-value, chromosome, position, reference and alternate alleles. They also contain complete variant annotation information, as well as results to the comparison of variants to GWAS Catalog.
The columns are decsribed in more detail in the following table:
Column name  |  Description
--- | ---
\#chrom  | chromosome of variant. from the original summary statistic.
pos	    | position coordinate of variant. from the original summary statistic.
ref	    | reference allele. from the original summary statistic.
alt	    | alternate allele. from the original summary statistic.
rsids	| variant RSID, if it exists. from the original summary statistic.
nearest_genes | from the original summary statistic.  
pval	| variant p-value. from the original summary statistic.
beta	| variant effect size. from the original summary statistic.
sebeta	| variant effect standard error. from the original summary statistic.
maf	    | variant minimum allele frequency. from the original summary statistic.
maf_case| variant minimum allele frequency in cases. from the original summary statistic.
maf_cont| variant minimum allele frequency in controls. from the original summary statistic.
variant| variant identifier.
locus_id| variant group identifier. Same as in the group reports.
pos_rmax| largest position in the variant group. Same as in group reports.
pos_rmin| smallest position in the variant group. Same as in group reports.
GENOME_* | Annotations from the gnoMAD genome information. Includes columns such as AF(allele frequency), FI_enrichment_nfe (Finnish enrichment agains Non-Finnish Europeans),FI_enrichment_nfe_est(Finnish enrichment agains Non-Finnish Europeans without Estonians).
EXOME_* | Annotations from the gnoMAD exome information. Includes columns such as AF(allele frequency), FI_enrichment_nfe (Finnish enrichment agains Non-Finnish Europeans),FI_enrichment_nfe_est(Finnish enrichment agains Non-Finnish Europeans without Estonians),FI_enrichment_nfe_swe(Finnish enrichment agains Non-Finnish Europeans without Swedish),FI_enrichment_nfe_est_swe(Finnish enrichment agains Non-Finnish Europeans without Estonians&Swedish populations).
most_severe_gene | The gene linked to the most severe consequence for this variant.
most_severe_consequence | most severe consequence linked to this variant.
FG_INFO_* | FinnGen INFO metric for the variant. In addition to the mean info score FG_INFO, the batch-specific info scores are reported on their own columns.
pval_trait | The p-value of a matching GWAS Catalog association. Empty if no match.
variant_hit | The variant identifier for a matching association from GWAS Catalog. Empty if no match.
trait | matching trait EFO code. Empty if no match.
trait_name | matching trait plaintext description. Empty if no match.

### Annotated Variants (.annotate.out)

These files contain the filtered, grouped and annotated variants. All of the information is in the full reports.
The columns are decsribed in more detail in the following table:
Column name  |  Description
--- | ---
\#chrom  | chromosome of variant. from the original summary statistic.
pos	    | position coordinate of variant. from the original summary statistic.
ref	    | reference allele. from the original summary statistic.
alt	    | alternate allele. from the original summary statistic.
rsids	| variant RSID, if it exists. from the original summary statistic.
nearest_genes | from the original summary statistic.  
pval	| variant p-value. from the original summary statistic.
beta	| variant effect size. from the original summary statistic.
sebeta	| variant effect standard error. from the original summary statistic.
maf	    | variant minimum allele frequency. from the original summary statistic.
maf_case| variant minimum allele frequency in cases. from the original summary statistic.
maf_cont| variant minimum allele frequency in controls. from the original summary statistic.
variant| variant identifier.
locus_id| variant group identifier. Same as in the group reports.
pos_rmax| largest position in the variant group. Same as in group reports.
pos_rmin| smallest position in the variant group. Same as in group reports.
GENOME_* | Annotations from the gnoMAD genome information. Includes columns such as AF(allele frequency), FI_enrichment_nfe (Finnish enrichment agains Non-Finnish Europeans),FI_enrichment_nfe_est(Finnish enrichment agains Non-Finnish Europeans without Estonians).
EXOME_* | Annotations from the gnoMAD exome information. Includes columns such as AF(allele frequency), FI_enrichment_nfe (Finnish enrichment agains Non-Finnish Europeans),FI_enrichment_nfe_est(Finnish enrichment agains Non-Finnish Europeans without Estonians),FI_enrichment_nfe_swe(Finnish enrichment agains Non-Finnish Europeans without Swedish),FI_enrichment_nfe_est_swe(Finnish enrichment agains Non-Finnish Europeans without Estonians&Swedish populations).
most_severe_gene | The gene linked to the most severe consequence for this variant.
most_severe_consequence | most severe consequence linked to this variant.
FG_INFO_* | FinnGen INFO metric for the variant. In addition to the mean info score FG_INFO, the batch-specific info scores are reported on their own columns.

### Filtered Variant Reports (.fetch.out)

These files contain the filtered and grouped variants. All of the information is in the full reports.
The columns are decsribed in more detail in the following table:
Column name  |  Description
--- | ---
\#chrom  | chromosome of variant. from the original summary statistic.
pos	    | position coordinate of variant. from the original summary statistic.
ref	    | reference allele. from the original summary statistic.
alt	    | alternate allele. from the original summary statistic.
rsids	| variant RSID, if it exists. from the original summary statistic.
nearest_genes | from the original summary statistic.  
pval	| variant p-value. from the original summary statistic.
beta	| variant effect size. from the original summary statistic.
sebeta	| variant effect standard error. from the original summary statistic.
maf	    | variant minimum allele frequency. from the original summary statistic.
maf_case| variant minimum allele frequency in cases. from the original summary statistic.
maf_cont| variant minimum allele frequency in controls. from the original summary statistic.
#variant| variant identifier.
locus_id| variant group identifier. Same as in the group reports.
pos_rmax| largest position in the variant group. Same as in group reports.
pos_rmin| smallest position in the variant group. Same as in group reports.
