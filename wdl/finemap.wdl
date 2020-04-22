import "finemap_sub.wdl" as sub

task preprocess {
    String pheno
    String sumstats_pattern
    File sumstats = sub(sumstats_pattern,"\\{PHENO\\}",pheno)
    String zones
    String docker
    Int cpu
    Int mem
    Boolean scale_se_by_pval
    Boolean x_chromosome
    String rsid_col
    String chromosome_col
    String position_col
    String allele1_col
    String allele2_col
    String freq_col
    String beta_col
    String se_col
    String p_col
    # can be helpful if adding finemapping with relaxed threshold after more stringent has already ben run.
    # does not include regions with lead snp < this
    Float p_threshold
    Float r2_threshold
    Float? minimum_pval

    command {

        make_finemap_inputs.py \
            --sumstats ${sumstats} \
            --rsid-col "${rsid_col}" \
            --chromosome-col "${chromosome_col}" \
            --position-col ${position_col} \
            --allele1-col ${allele1_col} \
            --allele2-col ${allele2_col} \
            --freq-col ${freq_col} \
            --beta-col ${beta_col} \
            --se-col ${se_col} \
            --p-col ${p_col} \
            --r2-threshold ${r2_threshold} \
            --set-rsid \
            --grch38 \
            --exclude-MHC \
            --no-upload \
            --gsdir '' \
            --localdir '' \
            --input-samples '' \
            --input-incl-samples '' \
            -n 1 \
            --var-y 1 \
            --out ${pheno} \
            ${true='--scale-se-by-pval ' false=' ' scale_se_by_pval} \
            ${true='--x-chromosome' false=' ' x_chromosome} \
            --p-threshold ${p_threshold} \
            ${true='--min-p-threshold ' false='' defined(minimum_pval)}${minimum_pval}

            res=`cat ${pheno}_had_results`

            if [ "$res" == "False" ]; then
                touch ${pheno}".z"
                touch ${pheno}".lead_snps.txt"
                touch ${pheno}".lead_snps.txt"
                touch ${pheno}".bed"
            fi
    }

    output {

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
        noAddress: false
    }
}

workflow finemap {

    String zones
    String docker
    String sumstats_pattern
    File phenolistfile
    File phenotypes

    Array[String] phenos = read_lines(phenolistfile)

    scatter (pheno in phenos) {

        call preprocess {
            input: zones=zones, docker=docker, pheno=pheno, sumstats_pattern=sumstats_pattern
        }

        if( preprocess.had_results) {
            call sub.ldstore_finemap {
                input: zones=zones, docker=docker, pheno=pheno, zfiles=preprocess.zfiles, phenofile=phenotypes, pheno=pheno
            }
        }

    }
}
