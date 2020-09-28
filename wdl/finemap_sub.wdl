task ldstore {
    String pheno
    String bgen_pattern
    File zfile
    File sample
    File phenofile
    String incl = pheno + ".incl"
    String prefix = basename(zfile, ".z")
    String chrom = sub(sub(prefix, pheno + "\\.chr", ""), "\\.[0-9\\-]+$", "")
    String bgenbucket = sub(sub(bgen_pattern, "^gs://", ""), "/.+$", "")
    String mountpoint = "/cromwell_root/gcsfuse/" + bgenbucket
    String bgen = sub(sub(bgen_pattern, "\\{CHR\\}", chrom), "^gs://" + bgenbucket, mountpoint)
    String bgen_gs = sub(bgen_pattern, "\\{CHR\\}", chrom)
    File bgi = sub(bgen_pattern, "\\{CHR\\}", chrom) + ".bgi"
    String master = prefix + ".master"
    String n_samples_file = prefix + ".n_samples.txt"
    String zones
    String docker
    Int cpu
    Int mem
    Boolean enable_fuse

    command <<<
        #!/usr/bin/env bash
        # mount bgen bucket

        mkdir -p ${mountpoint}

        if [[ ${enable_fuse} == "true" ]]
        then
            gcsfuse --implicit-dirs ${bgenbucket} ${mountpoint}
        else
            bgen_dir=$(dirname ${bgen})
            mkdir -p $bgen_dir
            gsutil -q cp ${bgen_gs} $bgen_dir/
        fi

        catcmd="cat"
        if [[ ${phenofile} == *.gz ]] || [[ ${phenofile} == *.bgz ]]
        then
            catcmd="zcat"
        fi
        echo "Reading phenotype file with $catcmd"

        $catcmd ${phenofile} | awk -v ph=${pheno} '
        BEGIN {
            FS = "\t"
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                h[$i] = i
            }
            exists = (ph in h)
            if (!exists) {
                print "Phenotype:"ph" not found in the given phenotype file." > "/dev/stderr"
                exit 1
            }
        }
        NR > 1 && $h[ph] != "NA" {
            print $1
        }
        ' > ${incl}

        if [[ $? -ne 0 ]]
        then
            echo "Given phenotype ${pheno} not found in given phenotype file ${phenofile}"
            exit 1
        fi

        wc -l ${incl} | cut -f1 -d' ' > ${n_samples_file}
        awk -v n_samples=`cat ${n_samples_file}` '
        BEGIN {
            OFS = ";"
            print "z", "bgen", "bgi", "bdose", "bcor", "ld", "sample", "incl", "n_samples"
            print "${zfile}", "${bgen}", "${bgi}", "${prefix}.bdose", "${prefix}.bcor", "${prefix}.ld", "${sample}", "${incl}", n_samples
        }' > ${master}

        n_threads=`grep -c ^processor /proc/cpuinfo`
        ldstore --in-files ${master} --write-bcor --read-only-bgen --n-threads $n_threads
        ldstore --in-files ${master} --bcor-to-text
        bgzip -@ $n_threads ${prefix}.ld
        mv ${prefix}.ld.gz ${prefix}.ld.bgz

        echo "Finished"

        if [[ ${enable_fuse} == "true" ]]
        then
            fusermount -u ${mountpoint}
        fi

        exit 0
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
        disks: "local-disk 200 HDD"
        zones: "${zones}"
        preemptible: 2
        noAddress: false
    }
}

task finemap {
    Int n_samples
    Int n_causal_snps
    File zfile
    File bcor
    File phenofile
    String pheno
    String prefix = basename(zfile, ".z")
    String ld = prefix + '.ld'
    String master = prefix + ".master"
    String zones
    String docker
    Int cpu=8
    Int mem=52

    command <<<
        #!/usr/bin/env bash
        awk -v n_samples=${n_samples} '
        BEGIN {
            OFS = ";"
            print "z", "bcor", "snp", "config", "cred", "n_samples", "log"
            print "${zfile}", "${bcor}", "${prefix}.snp", "${prefix}.config", "${prefix}.cred", n_samples, "${prefix}.log"
        }' > ${master}

		catcmd="cat"
		if [[ ${phenofile} == *.gz ]] || [[ ${phenofile} == *.bgz ]]
		then
			catcmd="zcat"
		fi

        prior_std=$($catcmd ${phenofile} | awk -v ph=${pheno} '
        BEGIN {
            FS = "\t"
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                h[$i] = i
            }
            exists=ph in h
            if (!exists) {
                print "Phenotype:"ph" not found in the given phenotype file." > "/dev/stderr"
                err = 1
                exit 1
            }
            cases = 0
            controls = 0
        }
        NR > 1 {
            vals[$(h[ph])] += 1
            if ($h[ph] != "NA" && $h[ph] != 0 && $h[ph] != 1) {
                print "Phenotype:"ph" seems a quantitative trait. Setting prior_std = 0.05."> "/dev/stderr"
                print 0.05
                err = 1
                exit 0
            }
        }
        END {
            if (!err) {
                phi = vals["1"]/(vals["1"]+vals["0"])
                std = 0.05 * sqrt(phi*(1-phi))
                printf std
            }
        }')

        if [[ $? -ne 0 ]]
        then
            echo "Error occurred while getting prior std from case control counts:" $prior_std
            exit 1
        fi

        n_threads=`grep -c ^processor /proc/cpuinfo`
        finemap --sss \
            --in-files ${master} \
            --log \
            --n-causal-snps ${n_causal_snps} \
            --n-threads $n_threads \
            --prior-std $prior_std

        # Merge p column
        cp ${prefix}.snp ${prefix}.snp.temp
        awk '
        BEGIN {
            OFS = " "
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                col[$i] = i
            }
        }
        FNR == NR {
            p[$col["rsid"]] = $col["p"]
            next
        }
        FNR < NR && FNR == 1 {
            print $0, "p"
        }
        FNR < NR && FNR > 1 {
            print $0, p[$2]
        }
        ' ${zfile} ${prefix}.snp.temp > ${prefix}.snp

    >>>

    output {
        File snp = prefix + ".snp"
        File config = prefix + ".config"
        # File cred = prefix + ".cred"
        Array[File] cred_files = glob("*.cred*")
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
    File zfile
    File ld_bgz
    File phenofile
    String pheno
    String prefix = basename(zfile, ".z")
    String zones
    String docker
    Int cpu=8
    Int mem=256
    Float min_cs_corr

    command <<<
        #!/usr/bin/env bash

		catcmd="cat"
		if [[ ${phenofile} == *.gz ]] || [[ ${phenofile} == *.bgz ]]
		then
			catcmd="zcat"
		fi

        var_y=$($catcmd ${phenofile} | awk -v ph=${pheno} '
        BEGIN {
            FS = "\t"
        }
        NR == 1 {
            for(i = 1; i <= NF; i++) {
                h[$i] = i
            }
            exists=ph in h
            if (!exists) {
                print "Phenotype:"ph" not found in the given phenotype file." > "/dev/stderr"
                err = 1
                exit 1
            }
            cases = 0
            controls = 0
        }
        NR > 1 {
            vals[$h[ph]] += 1
            if ($h[ph] != "NA" && $h[ph] != 0 && $h[ph] != 1) {
                print "Phenotype:"ph" seems a quantitative trait. Setting var_y = 1." > "/dev/stderr"
                print 1.0
                err = 1
                exit 0
            }
        }
        END {
            if (!err) {
                phi = vals["1"] / (vals["1"]+vals["0"])
                var_y = phi*(1-phi)
                printf var_y
            }
        }')

        if [[ $? -ne 0 ]]
        then
            echo "Error occurred while getting var_y from case control counts:" $var_y
            exit 1
        fi

        run_susieR.R \
            --z ${zfile} \
            --ld ${ld_bgz} \
            -n ${n_samples} \
            --L ${n_causal_snps} \
            --var-y $var_y \
            --snp ${prefix}.susie.snp \
            --cred ${prefix}.susie.cred \
            --log ${prefix}.susie.log \
            --susie-obj ${prefix}.susie.rds \
            --save-susie-obj \
            --write-alpha \
            --write-single-effect \
            --min-cs-corr ${min_cs_corr}

    >>>

    output {
        File log = prefix + ".susie.log"
        File snp = prefix + ".susie.snp"
        File cred = prefix + ".susie.cred"
        File rds = prefix + ".susie.rds"

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
                col[$i] = i
            }
            gsub(" ", "\t")
            print "trait", "region", "v", $0
        }
        FNR == 1 {
            match(FILENAME, /(chr[0-9X]+)\.([0-9]+-[0-9]+)\./, a)
            region = a[1]":"a[2]
        }
        FNR > 1 {
            chrom = substr($col["chromosome"], 4)
            sub(/^0/, "", chrom)
            v = sprintf( \
                "%s:%s:%s:%s", \
                chrom, \
                $col["position"], \
                $col["allele1"], \
                $col["allele2"] \
            )
            gsub(" ", "\t")
            print pheno, region, v, $0 | "sort -V -k2,3"
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
                col[$i] = i
            }
            gsub(" ", "\t")
            print "trait", "region", $0
        }
        FNR == 1 {
            match(FILENAME, /(chr[0-9X]+)\.([0-9]+-[0-9]+)\./, a)
            region = a[1]":"a[2]
        }
        FNR > 1 {
            gsub(" ", "\t")
            print pheno, region, $0 | "sort -V -k2,3"
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
                if ($i != "NA") {
                    print $i, cred, i/2
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
            a[$1":"$2] = $3
            next
        }
        FNR < NR && FNR == 1 {
            for (i = 1; i <= NF; i++) {
                col[$i] = i
            }
            cs_str = "cs"
            for (i = 2; i <= ${n_causal_snps}; i++) {
                cs_str = cs_str"\tcs"i
            }
            print $0, cs_str
        }
        FNR < NR && FNR > 1 {
            cs_str = ""
            for (i = 1; i <= ${n_causal_snps}; i++) {
                cs_str = cs_str"\t"(($col["rsid"]":"i in a) ? a[$col["rsid"]":"i] : "-1")
            }
            print $0""cs_str
        }
        ' ${pheno}.FINEMAP.temp.cred ${pheno}.FINEMAP.temp.snp | bgzip -c -@ ${cpu} > ${pheno}.FINEMAP.snp.bgz
        tabix -s 6 -b 7 -e 7 -S 1 ${pheno}.FINEMAP.snp.bgz

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
            match(FILENAME, /(chr[0-9X]+)\.([0-9]+-[0-9]+)\./, a)
            region = a[1]":"a[2]
            regions[region] = region
        }
        FNR > 1 && $0 ~ /Regional SNP heritability/ {
            match($0, /([0-9\.\-e]+) \(SD: ([0-9\.\-e]+) ; 95% CI: \[([0-9\.\-e]+),([0-9\.\-e]+)/, a)
            h2g[region] = a[1]"\t"a[2]"\t"a[3]"\t"a[4]
        }
        FNR > 1 && $0 ~ /Log10-BF of >= one causal SNP/ {
            match($0, /: ([0-9\.\-]+)/, a)
            log10bf[region] = a[1]
        }
        FNR > 1 && $0 ~ /Post-expected # of causal SNPs/ {
            match($0, /: ([0-9\.\-]+)/, a)
            n_exp[region] = a[1]
        }
        FNR > 1 && $0 ~ /Post-Pr\(# of causal SNPs is k\)/ {
            getline
            pp_str=""
            for (i = 1; i <= ${n_causal_snps}; i++) {
                getline
                match($0, /-> ([0-9\.]+)/, a)
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
        tabix -s 5 -b 6 -e 6 -S 1 ${pheno}.SUSIE.snp.bgz

        # Combine susie .cred files
        awk -v pheno=${pheno} '
        BEGIN {
            OFS = "\t"
        }
        NR == 1 {
            print "trait", "region", $0
        }
        FNR == 1 {
            match(FILENAME, /(chr[0-9X]+)\.([0-9]+-[0-9]+)\./, a)
            region = a[1]":"a[2]
        }
        FNR > 1 {
            print pheno, region, $0
        }
        ' ${sep=" " susie_cred} | bgzip -c -@ ${cpu} > ${pheno}.SUSIE.cred.bgz

    >>>

    output {
        File out_finemap_snp = pheno + ".FINEMAP.snp.bgz"
        File out_finemap_snp_tbi = pheno + ".FINEMAP.snp.bgz.tbi"
        File out_finemap_config = pheno + ".FINEMAP.config.bgz"
        File out_finemap_region = pheno + ".FINEMAP.region.bgz"
        File out_susie_snp = pheno + ".SUSIE.snp.bgz"
        File out_susie_snp_tbi = pheno + ".SUSIE.snp.bgz.tbi"
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


task filter_and_summarize{
    String zones
    String docker
    Int cpu
    Int mem
    String pheno
    File susie_snps
    File susie_snps_tbi
    File susie_cred
    File? snp_annot_file
    File? snp_annot_file_tbi
    String? snp_annotation_fields
    String? set_variant_id_map_chr
    Float good_cred_r2
    Int mock
    command <<<
        filter_and_summarize.py --min_r2 ${good_cred_r2} ${susie_cred} \
            ${susie_snps} ${pheno}.SUSIE \
            ${true='--variant_annot ' false='' defined(snp_annot_file)}${snp_annot_file} \
            ${true='--variant_annot_cols ' false='' defined(snp_annotation_fields)}${snp_annotation_fields} \
            ${true='--set_variant_id_map_chr ' false='' defined(set_variant_id_map_chr)}${set_variant_id_map_chr}
    >>>

    output {
        File out_susie_snp_filtered = pheno + ".SUSIE.snp.filter.tsv"
        File out_susie_cred_summary = pheno + ".SUSIE.cred.summary.tsv"
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
    File phenofile

    String? set_variant_id_map_chr

    scatter (zfile in zfiles) {

        call ldstore {
            input: zones=zones, docker=docker, pheno=pheno, phenofile=phenofile, zfile=zfile
        }

        call finemap {
            input: zones=zones, docker=docker, zfile=zfile,
                bcor=ldstore.bcor, n_samples=ldstore.n_samples,
                n_causal_snps=n_causal_snps, phenofile=phenofile, pheno=pheno
        }

        call susie {
            input: zones=zones, zfile=zfile, ld_bgz=ldstore.ld_bgz, n_samples=ldstore.n_samples, n_causal_snps=n_causal_snps,
                phenofile=phenofile, pheno=pheno
        }
    }

    call combine {
        input: zones=zones, docker=docker, pheno=pheno, n_causal_snps=n_causal_snps,
            finemap_config=finemap.config, finemap_snp=finemap.snp, finemap_cred_files=finemap.cred_files, finemap_log=finemap.log,
            susie_snp=susie.snp, susie_cred=susie.cred
    }

    call filter_and_summarize {
        input: zones=zones, pheno=pheno, docker=docker, susie_snps=combine.out_susie_snp,
            susie_snps_tbi=combine.out_susie_snp_tbi, susie_cred=combine.out_susie_cred,
            set_variant_id_map_chr=set_variant_id_map_chr, mock=1
    }

}
