import "finemap_tasks.wdl" as tasks


workflow finemap {

    String zones
    String docker
    File zfile_list_file
    # File phenotype_summary_file
    Int n_causal_snps

    Array[String] zfiles = read_lines(zfile_list_file)
    # Array[Array[String]] phenotype_summary = read_tsv(phenotype_summary_file)
    # String pheno = phenotype_summary[0]
    String pheno = "E4_HYTHY_AI_STRICT_EXTRA"
    Int n_samples = 397278
    Int n_cases = 26966
    Float phi = 1.0 * n_cases / n_samples
    Float var_y = phi * (1 - phi)

    scatter (zfile in zfiles) {

        File ld_bgz = sub(zfile, ".z$", ".ld.bgz")

        call tasks.finemap as finemap {
            input: zones=zones, docker=docker, zfile=zfile,
                ld_bgz=ld_bgz, n_samples=n_samples,
                n_causal_snps=n_causal_snps, var_y = var_y, pheno=pheno
        }

        call tasks.susie as susie {
            input: zones=zones, zfile=zfile, ld_bgz=ld_bgz, n_samples=n_samples, n_causal_snps=n_causal_snps, var_y = var_y,
                pheno=pheno
        }
    }

    call tasks.combine as combine {
        input: zones=zones, docker=docker, pheno=pheno, n_causal_snps=n_causal_snps,
            finemap_config=finemap.config, finemap_snp=finemap.snp, finemap_cred_files=finemap.cred_files, finemap_log=finemap.log,
            susie_snp=susie.snp, susie_cred=susie.cred
    }

    call tasks.filter_and_summarize as filter_and_summarize {
        input: zones=zones, pheno=pheno, docker=docker, susie_snps=combine.out_susie_snp,
            susie_snps_tbi=combine.out_susie_snp_tbi, susie_cred=combine.out_susie_cred
    }
}
