import "finemap_tasks.wdl" as tasks

task preprocess {
    String pheno
    String? prefix
    String sumstats_pattern
    File sumstats = sub(sumstats_pattern,"\\{PHENO\\}",pheno)
    String zones
    String docker
    Int cpu
    Int mem
    Boolean scale_se_by_pval
    Boolean x_chromosome
    Boolean set_variant_id
    String rsid_col
    String chromosome_col
    String position_col
    String allele1_col
    String allele2_col
    String freq_col
    String beta_col
    String se_col
    String p_col
    String delimiter
    # can be helpful if adding finemapping with relaxed threshold after more stringent has already ben run.
    # does not include regions with lead snp < this
    Float p_threshold
    Float? minimum_pval
    Int? window

    command <<<

        zcat "${sumstats}" | awk '
        BEGIN {
            OFS = "\t"
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                col[$i] = i
            }
            chr = "${chromosome_col}"
            pos = "${position_col}"
            ref = "${allele1_col}"
            alt = "${allele2_col}"
            beta = "${beta_col}"
            se = "${se_col}"
            pval = "${p_col}"
            print chr, pos, ref, alt, "maf", beta, se, pval
        }
        NR > 1 {
            if ($col[chr] == 23) {
                $col[chr] = "X"
            }
            # avoid duplication for now
            if ($col["SNP"] in snp) {
                next
            } else {
                snp[$col["SNP"]] = 1
            }
            af = ($col["E4_HYTHY_AI_STRICT_C6_AF_Allele2"] * $col["E4_HYTHY_AI_STRICT_C6_N"] + $col["UKBB_EUR_E4_HYTHY_AI_STRICT_EXTRA_AF_Allele2"] * $col["UKBB_EUR_E4_HYTHY_AI_STRICT_EXTRA_N"] + $col["UKBB_AFR_E4_HYTHY_AI_STRICT_EXTRA_AF_Allele2"] * $col["UKBB_AFR_E4_HYTHY_AI_STRICT_EXTRA_N"] + $col["UKBB_CSA_E4_HYTHY_AI_STRICT_EXTRA_AF_Allele2"] * $col["UKBB_CSA_E4_HYTHY_AI_STRICT_EXTRA_N"]) / ($col["E4_HYTHY_AI_STRICT_C6_N"] + $col["UKBB_EUR_E4_HYTHY_AI_STRICT_EXTRA_N"] + $col["UKBB_AFR_E4_HYTHY_AI_STRICT_EXTRA_N"] + $col["UKBB_CSA_E4_HYTHY_AI_STRICT_EXTRA_N"])
            print $col[chr], $col[pos], $col[ref], $col[alt], af, $col[beta], $col[se], $col[pval]
        }
        ' > "${pheno}.tsv"

        head "${pheno}.tsv"

        make_finemap_inputs.py \
            --sumstats "${pheno}.tsv" \
            --rsid-col "${rsid_col}" \
            --chromosome-col "${chromosome_col}" \
            --position-col "${position_col}" \
            --allele1-col "${allele1_col}" \
            --allele2-col "${allele2_col}" \
            --freq-col "${freq_col}" \
            --beta-col "${beta_col}" \
            --se-col "${se_col}" \
            --p-col "${p_col}" \
            --delimiter "${delimiter}" \
            --grch38 \
            --no-upload \
            --prefix ${pheno} \
            --out ${pheno} \
            --wdl \
            ${true='--scale-se-by-pval ' false=' ' scale_se_by_pval} \
            ${true='--x-chromosome' false=' ' x_chromosome} \
            ${true='--set-variant-id ' false=' ' set_variant_id} \
            ${true='--variant-id-sep ":" ' false=' ' set_variant_id} \
            --p-threshold ${p_threshold} \
            ${true='--min-p-threshold ' false='' defined(minimum_pval)}${minimum_pval} \
            ${true="--window " false=' ' defined(window)}${window}

            res=`cat ${pheno}_had_results`

            if [ "$res" == "False" ]; then
                touch ${pheno}".z"
                touch ${pheno}".lead_snps.txt"
                touch ${pheno}".lead_snps.txt"
                touch ${pheno}".bed"
            fi
    >>>

    output {

        File formatted_tsv = pheno + ".tsv"
        Array[File] zfiles = glob("*.z")
        File leadsnps = pheno + ".lead_snps.txt"
        File bed = pheno + ".bed"
        File log = pheno + ".log"
        Boolean had_results = read_boolean("${pheno}_had_results")
    }

    runtime {
        docker: "${docker}"
        cpu: "${cpu}"
        memory: "${mem} GB"
        disks: "local-disk 20 HDD"
        zones: "${zones}"
        preemptible: 2
        noAddress: true
    }
}

task postfilter {
    Array[File] in_zfiles
    String zones
    String docker
    Int cpu = 1
    Int mem = 1
    command <<<
        for file in ${sep=" " in_zfiles}
        do
            if [[ $(wc -l $file | cut -f1 -d' ') -lt 64142 ]]
            then
                cp $file $(basename $file)
            fi

        done
    >>>

    output {

        Array[File] zfiles = glob("*.z")
    }

    runtime {
        docker: "${docker}"
        cpu: "${cpu}"
        memory: "${mem} GB"
        disks: "local-disk 20 HDD"
        zones: "${zones}"
        preemptible: 2
        noAddress: true
    }
}


task export_ld {
    Array[String] pops
    Array[Int] n_samples
    File zfile_meta
    String prefix = basename(zfile_meta, ".z")
    String zones
    String docker
    Int cpu
    Int mem

    File script

    command <<<

        PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory=${mem}g pyspark-shell" \
        TMPDIR=/cromwell_root/ \
        python3 ${script} \
        --pops ${sep=" " pops} \
        --n-samples ${sep=" " n_samples} \
        --zfile-meta "${zfile_meta}" \
        --meta \
        --sparsity full \
        --tmpdir $TMPDIR \
        --out ${prefix} \
        && touch _SUCCESS

    >>>

    output {

        File success = "_SUCCESS"
        File ld_bgz = prefix + ".ld.bgz"
        File zfile = prefix + ".z"

    }

    runtime {

        docker: "${docker}"
        cpu: "${cpu}"
        memory: "102 GB"
        disks: "local-disk 200 HDD"
        zones: "${zones}"
        preemptible: 2
        noAddress: false
    }
}

workflow finemap {

    String zones
    String docker
    # File zfile_list_file
    # File phenotype_summary_file
    Array[String] pops

    # Array[String] zfiles = read_lines(zfile_list_file)
    # Array[Array[String]] phenotype_summary = read_tsv(phenotype_summary_file)
    # String pheno = phenotype_summary[0]
    String sumstats_pattern
    String pheno
    Array[Int] n_samples
    Map[String, Int] n_samples_dict
    Int n_causal_snps
    Float var_y

    call preprocess {
        input: zones=zones, docker=docker, pheno=pheno, sumstats_pattern=sumstats_pattern
    }

    call postfilter {
        input: zones=zones, docker=docker, in_zfiles=preprocess.zfiles
    }

    # File zfile = preprocess.zfiles[12]
    scatter (zfile in postfilter.zfiles) {

        call export_ld {
            input: zones=zones, pops=pops, n_samples=n_samples, zfile_meta=zfile
        }

        call tasks.finemap as finemap {
            input: zones=zones, docker=docker, zfile=export_ld.zfile,
                ld_bgz=export_ld.ld_bgz, n_samples=n_samples_dict["meta"],
                n_causal_snps=n_causal_snps, pheno=pheno, var_y=var_y
        }

        call tasks.susie as susie {
            input: zones=zones, zfile=export_ld.zfile, ld_bgz=export_ld.ld_bgz, n_samples=n_samples_dict["meta"], n_causal_snps=n_causal_snps,
                pheno=pheno, var_y=var_y
        }

    }

    call tasks.combine as combine {
        input: zones=zones, docker=docker, pheno=pheno, n_causal_snps=n_causal_snps,
            finemap_config=finemap.config, finemap_snp=finemap.snp, finemap_cred_files=finemap.cred_files, finemap_log=finemap.log,
            susie_snp=susie.snp, susie_cred=susie.cred
    }

}
