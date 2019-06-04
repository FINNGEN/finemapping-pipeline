import "finemap_sub.wdl" as sub

task preprocess {
    String pheno
    String sumstatsdir
    File sumstats = sumstatsdir + "/" + pheno + ".gz"
    String zones
    String docker
    Int cpu
    Int mem
    Boolean scale_se_by_pval

    command {

        make_finemap_inputs.py \
            --sumstats ${sumstats} \
            --rsid-col "SNPID" \
            --chromosome-col "CHR" \
            --position-col "POS" \
            --allele1-col "Allele1" \
            --allele2-col "Allele2" \
            --freq-col "AF_Allele2" \
            --beta-col "BETA" \
            --se-col "SE" \
            --p-col "p.value" \
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
            ${true='--scale-se-by-pval' false=' ' scale_se_by_pval}
    }

    output {

        Array[File] zfiles = glob("*.z")
        File leadsnps = pheno + ".lead_snps.txt"
        File bed = pheno + ".bed"
        File log = pheno + ".log"

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

workflow finemap {

    String zones
    String docker
    String sumstatsdir
    File phenolistfile
    File phenotypes

    Array[String] phenos = read_lines(phenolistfile)

    scatter (pheno in phenos) {

        call preprocess {
            input: zones=zones, docker=docker, pheno=pheno, sumstatsdir=sumstatsdir
        }

        call sub.ldstore_finemap {
            input: zones=zones, docker=docker, pheno=pheno, zfiles=preprocess.zfiles, phenofile=phenotypes, pheno=pheno
        }
    }
}
