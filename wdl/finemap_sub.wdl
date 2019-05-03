task ldstore {
    String pheno
    String incldir
    String bgendir
    File zfile
    File sample
    File incl = incldir + "/" + pheno + ".incl"
    String prefix = basename(zfile, ".z")
    String chrom = sub(sub(prefix, pheno + "\\.chr", ""), "\\.[0-9\\-]+$", "")
    File bgen = bgendir + "/" + chrom + "R3.bgen"
    File bgi = bgen + ".bgi"
    String master = prefix + ".master"
    String n_samples_file = prefix + ".n_samples.txt"
    String docker
    Int cpu
    Int mem

    command <<<

        wc -l ${incl} | cut -f1 -d' ' > ${n_samples_file}
        cat ${n_samples_file}
        awk -v n_samples=`cat ${n_samples_file}` '
        BEGIN {
            OFS = ";"
            print "z", "bgen", "bgi", "bcor", "ld", "sample", "incl", "n_samples"
            print "${zfile}", "${bgen}", "${bgi}", "${prefix}.bcor", "${prefix}.ld", "${sample}", "${incl}", n_samples
        }' > ${master}

        cat ${master}

        ldstore --in-files ${master} --write-bcor --n-threads ${cpu}
        ldstore --in-files ${master} --bcor-to-text

        bgzip -@ ${cpu} ${prefix}.ld
        mv ${prefix}.ld.gz ${prefix}.ld.bgz

    >>>

    output {

        Int n_samples = read_int(n_samples_file)
        File bcor = prefix + ".bcor"
        File ld_bgz = prefix + ".ld.bgz"

    }

    runtime {

        docker: "${docker}"
        cpu: "${cpu}"
        memory: "${mem} GB"
        disks: "local-disk 100 HDD"
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: false
    }
}

task finemap {
    Int n_samples
    Int n_causal_snps
    Float corr_group
    File zfile
    # File bcor
    File ld_bgz
    String prefix = basename(zfile, ".z")
    String master = prefix + ".master"
    String ld = basename(ld_bgz, ".bgz")
    String docker
    Int cpu
    Int mem

    command <<<

        zcat ${ld_bgz} > ${ld}
        awk -v n_samples=${n_samples} '
        BEGIN {
            OFS = ";"
            print "z", "ld", "snp", "config", "cred", "n_samples", "log"
            print "${zfile}", "${ld}", "${prefix}.snp", "${prefix}.config", "${prefix}.cred", n_samples, "${prefix}.log"
        }' > ${master}

        finemap --sss \
            --in-files ${master} \
            --log \
            --n-causal-snps ${n_causal_snps} \
            --group-snps \
            --corr-group ${corr_group}

    >>>

    output {

        File snp = prefix + ".snp"
        File config = prefix + ".config"
        File cred = prefix + ".cred"
        File log = prefix + ".log_sss"

    }

    runtime {

        docker: "${docker}"
        cpu: "${cpu}"
        memory: "${mem} GB"
        disks: "local-disk 25 HDD"
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: true
    }
}


task susie {
    Int n_samples
    Int n_causal_snps
    Float var_y
    File zfile
    File ld_bgz
    String prefix = basename(zfile, ".z")
    String docker
    Int cpu
    Int mem

    command {

        run_susieR.R \
            --z ${zfile} \
            --ld ${ld_bgz} \
            -n ${n_samples} \
            --L ${n_causal_snps} \
            --var-y ${var_y} \
            --snp ${prefix}.susie.snp \
            --cred ${prefix}.susie.cred \
            --log ${prefix}.susie.log

    }

    output {

        File snp = prefix + ".susie.snp"
        File cred = prefix + ".susie.cred"
        File log = prefix + ".susie.log"

    }

    runtime {

        docker: "${docker}"
        cpu: "${cpu}"
        memory: "${mem} GB"
        disks: "local-disk 50 HDD"
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: true
    }
}

# task combine {
#     String pheno
#     Array[File] finemap_config
#     Array[File] finemap_snp
#     Array[File] finemap_cred
#     Array[File] finemap_log
#     Array[File] susie_snp
#     Array[File] susie_cred
#     String docker
#     Int cpu
#     Int mem

#     command {

#         Rscript run_susieR.R \
#             --z ${zfile} \
#             --ld ${ld_bgz} \
#             -n ${n_samples} \
#             -L ${n_causal_snps} \
#             --var-y ${var_y} \
#             --snp ${prefix} + ".susie.snp" \
#             --cred ${prefix} + ".susie.cred" \
#             --log ${prefix} + ".susie.log"

#     }

#     output {

#         File snp = prefix + ".susie.snp"
#         File cred = prefix + ".susie.cred"
#         File log = prefix + ".susie.log"

#     }

#     runtime {

#         docker: "${docker}"
#         cpu: "${cpu}"
#         memory: "${mem} GB"
#         disks: "local-disk 25 HDD"
#         zones: "europe-west1-b"
#         preemptible: 2
#         noAddress: true
#     }
# }

workflow ldstore_finemap {

    String docker
    String pheno
    Int n_causal_snps
    Array[File] zfiles

    scatter (zfile in zfiles) {

        call ldstore {
            input: docker=docker, pheno=pheno, zfile=zfile
        }

        call finemap {
            input: docker=docker, zfile=zfile, ld_bgz=ldstore.ld_bgz, n_samples=ldstore.n_samples, n_causal_snps=n_causal_snps
        }

        call susie {
            input: zfile=zfile, ld_bgz=ldstore.ld_bgz, n_samples=ldstore.n_samples, n_causal_snps=n_causal_snps
        }
    }

    # call combine {
    #     input: docker=docker, pheno=pheno,
    #         finemap_config=finemap.config, finemap_snp=finemap.snp, finemap_cred=finemap.cred, finemap_log=finemap.log,
    #         susie_snp=susie.snp, susie_cred=susie.cred
    # }
}