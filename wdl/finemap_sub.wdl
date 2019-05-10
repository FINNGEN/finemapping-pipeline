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
    String zones
    String docker
    Int cpu
    Int mem

    command <<<

        wc -l ${incl} | cut -f1 -d' ' > ${n_samples_file}
        awk -v n_samples=`cat ${n_samples_file}` '
        BEGIN {
            OFS = ";"
            print "z", "bgen", "bgi", "bcor", "ld", "sample", "incl", "n_samples"
            print "${zfile}", "${bgen}", "${bgi}", "${prefix}.bcor", "${prefix}.ld", "${sample}", "${incl}", n_samples
        }' > ${master}

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
        zones: "${zones}"
        preemptible: 2
        noAddress: false
    }
}

task finemap {
    Int n_samples
    Int n_causal_snps
    Float corr_group
    File zfile
    File bcor
    String prefix = basename(zfile, ".z")
    String ld = prefix + '.ld'
    String master = prefix + ".master"
    String dollar = "$"
    String zones
    String docker
    Int cpu
    Int mem

    command <<<

        awk -v n_samples=${n_samples} '
        BEGIN {
            OFS = ";"
            print "z", "bcor", "snp", "config", "cred", "n_samples", "log"
            print "${zfile}", "${bcor}", "${prefix}.snp", "${prefix}.config", "${prefix}.cred", n_samples, "${prefix}.log"
        }' > ${master}

        finemap --sss \
            --in-files ${master} \
            --log \
            --n-causal-snps ${n_causal_snps} \
            --n-threads ${cpu}

        # Merge p column
        cp ${prefix}.snp ${prefix}.snp.temp
        awk '
        BEGIN {
            OFS = " "
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                col[${dollar}i] = i
            }
        }
        FNR == NR {
            p[${dollar}col["rsid"]] = ${dollar}col["p"]
            next
        }
        FNR < NR && FNR == 1 {
            print ${dollar}0, "p"
        }
        FNR < NR && FNR > 1 {
            print ${dollar}0, p[${dollar}2]
        }
        ' ${zfile} ${prefix}.snp.temp > ${prefix}.snp

    >>>

    output {

        File snp = prefix + ".snp"
        File config = prefix + ".config"
        # File cred = prefix + ".cred"
        Array[File] cred_files = glob(prefix + ".cred*")
        File log = prefix + ".log_sss"

    }

    runtime {

        docker: "${docker}"
        cpu: "${cpu}"
        memory: "${mem} GB"
        disks: "local-disk 100 HDD"
        zones: "${zones}"
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
    String zones
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
        zones: "${zones}"
        preemptible: 2
        noAddress: true
    }
}

task combine {
    String pheno
    Int n_causal_snps
    Array[File] finemap_config
    Array[File] finemap_snp
    Array[Array[File]] finemap_cred_files
    Array[File] finemap_cred = flatten(finemap_cred_files)
    Array[File] finemap_log
    Array[File] susie_snp
    Array[File] susie_cred
    String dollar = "$"
    String zones
    String docker
    Int cpu
    Int mem

    command <<<

        cat << "__EOF__" > combine_snp.awk
        BEGIN {
            OFS = "\t"
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                col[${dollar}i] = i
            }
            gsub(" ", "\t")
            print "trait", "region", "v", ${dollar}0
        }
        FNR == 1 {
            match(FILENAME, /(chr[0-9]+)\.([0-9]+-[0-9]+)\./, a)
            region = a[1]":"a[2]
        }
        FNR > 1 {
            v = sprintf( \
                "%s:%s:%s:%s", \
                int(substr(${dollar}col["chromosome"], 4)), \
                ${dollar}col["position"], \
                ${dollar}col["allele1"], \
                ${dollar}col["allele2"] \
            )
            gsub(" ", "\t")
            print pheno, region, v, ${dollar}0
        }
        __EOF__

        # Combine finemap .snp files
        awk -f combine_snp.awk -v pheno=${pheno} ${sep=" " finemap_snp} > ${pheno}.FINEMAP.temp.snp

        # Combine finemap .config files
        awk -v pheno=${pheno} '
        BEGIN {
            OFS = "\t"
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                col[${dollar}i] = i
            }
            gsub(" ", "\t")
            print "trait", "region", ${dollar}0
        }
        FNR == 1 {
            match(FILENAME, /(chr[0-9]+)\.([0-9]+-[0-9]+)\./, a)
            region = a[1]":"a[2]
        }
        FNR > 1 {
            gsub(" ", "\t")
            print pheno, region, ${dollar}0
        }
        ' ${sep=" " finemap_config} | bgzip -c -@ ${cpu} > ${pheno}.FINEMAP.config.bgz

        # Extract SNPs in finemap .cred files
        awk '
        BEGIN {
            OFS = "\t"
        }
        FNR == 2 {
            n_cols = NF
            match(FILENAME, "\\.cred([0-9]+)", a)
            cred = a[1]
        }
        FNR >= 4 {
            for (i = 2; i < n_cols; i += 2) {
                if (${dollar}i != "NA") {
                    print ${dollar}i, cred, i/2
                }
            }
        }
        ' ${sep=" " finemap_cred} > ${pheno}.FINEMAP.temp.cred

        # Merge with finemap .snp
        awk '
        BEGIN {
            OFS = "\t"
        }
        FNR == NR {
            a[${dollar}1":"${dollar}2] = ${dollar}3
            next
        }
        FNR < NR && FNR == 1 {
            for (i = 1; i <= NF; i++) {
                col[${dollar}i] = i
            }
            cs_str = "cs"
            for (i = 2; i <= ${n_causal_snps}; i++) {
                cs_str = cs_str"\tcs"i
            }
            print ${dollar}0, cs_str
        }
        FNR < NR && FNR > 1 {
            cs_str = ""
            for (i = 1; i <= ${n_causal_snps}; i++) {
                cs_str = cs_str"\t"((${dollar}col["rsid"]":"i in a) ? a[${dollar}col["rsid"]":"i] : "-1")
            }
            print ${dollar}0""cs_str
        }
        ' ${pheno}.FINEMAP.temp.cred ${pheno}.FINEMAP.temp.snp | bgzip -c -@ ${cpu} > ${pheno}.FINEMAP.snp.bgz

        # Extract region statistics from finemap .log_sss files
        awk -v pheno=${pheno} '
        BEGIN {
            OFS = "\t"
            pp_str = "prob_1SNP"
            for (i = 2; i <= ${n_causal_snps}; i++) {
                pp_str = pp_str"\tprob_"i"SNP"
            }
            print "trait", "region", "h2g", "h2g_sd", "h2g_lower95", "h2g_upper95", "log10bf", pp_str, "expectedvalue"
        }
        FNR == 1 {
            match(FILENAME, /(chr[0-9]+)\.([0-9]+-[0-9]+)\./, a)
            region = a[1]":"a[2]
            regions[region] = region
        }
        FNR > 1 && ${dollar}0 ~ /Regional SNP heritability/ {
            match(${dollar}0, /([0-9\.]+) \(SD: ([0-9\.]+) ; 95% CI: \[([0-9\.]+),([0-9\.]+)/, a)
            h2g[region] = a[1]"\t"a[2]"\t"a[3]"\t"a[4]
        }
        FNR > 1 && ${dollar}0 ~ /Log10-BF of >= one causal SNP/ {
            match(${dollar}0, /: ([0-9\.\-]+)/, a)
            log10bf[region] = a[1]
        }
        FNR > 1 && ${dollar}0 ~ /Post-expected # of causal SNPs/ {
            match(${dollar}0, /: ([0-9\.\-]+)/, a)
            n_exp[region] = a[1]
        }
        FNR > 1 && ${dollar}0 ~ /Post-Pr\(# of causal SNPs is k\)/ {
            getline
            pp_str=""
            for (i = 1; i <= ${n_causal_snps}; i++) {
                getline
                match(${dollar}0, /-> ([0-9\.]+)/, a)
                pp_str = pp_str""a[1]"\t"
            }
            pp[region] = pp_str""n_exp[region]
        }
        END {
            for (region in regions) {
                print pheno, region, h2g[region], log10bf[region], pp[region]
            }
        }
        ' ${sep=" " finemap_log} | bgzip -c -@ ${cpu} > ${pheno}.FINEMAP.region.bgz

        # Combine susie .snp files
        awk -f combine_snp.awk -v pheno=${pheno} ${sep=" " susie_snp} | bgzip -c -@ ${cpu} > ${pheno}.SUSIE.snp.bgz

        # Combine susie .cred files
        awk -v pheno=${pheno} '
        BEGIN {
            OFS = "\t"
        }
        NR == 1 {
            print "trait", "region", ${dollar}0
        }
        FNR == 1 {
            match(FILENAME, /(chr[0-9]+)\.([0-9]+-[0-9]+)\./, a)
            region = a[1]":"a[2]
        }
        FNR > 1 {
            print pheno, region, ${dollar}0
        }
        ' ${sep=" " susie_cred} | bgzip -c -@ ${cpu} > ${pheno}.SUSIE.cred.bgz

    >>>

    output {

        File out_finemap_snp = pheno + ".FINEMAP.snp.bgz"
        File out_finemap_config = pheno + ".FINEMAP.config.bgz"
        File out_finemap_region = pheno + ".FINEMAP.region.bgz"
        File out_susie_snp = pheno + ".SUSIE.snp.bgz"
        File out_susie_cred = pheno + ".SUSIE.cred.bgz"

    }

    runtime {

        docker: "${docker}"
        cpu: "${cpu}"
        memory: "${mem} GB"
        disks: "local-disk 30 HDD"
        zones: "${zones}"
        preemptible: 2
        noAddress: true
    }
}

workflow ldstore_finemap {

    String zones
    String docker
    String pheno
    Int n_causal_snps
    Array[File] zfiles

    scatter (zfile in zfiles) {

        call ldstore {
            input: zones=zones, docker=docker, pheno=pheno, zfile=zfile
        }

        call finemap {
            input: zones=zones, docker=docker, zfile=zfile, bcor=ldstore.bcor, n_samples=ldstore.n_samples, n_causal_snps=n_causal_snps
        }

        call susie {
            input: zones=zones, zfile=zfile, ld_bgz=ldstore.ld_bgz, n_samples=ldstore.n_samples, n_causal_snps=n_causal_snps
        }
    }

    call combine {
        input: zones=zones, docker=docker, pheno=pheno, n_causal_snps=n_causal_snps,
            finemap_config=finemap.config, finemap_snp=finemap.snp, finemap_cred_files=finemap.cred_files, finemap_log=finemap.log,
            susie_snp=susie.snp, susie_cred=susie.cred
    }
}