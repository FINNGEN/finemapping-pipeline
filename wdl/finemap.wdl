import "finemap_sub.wdl" as sub

task preprocess {
    String pheno
    File phenofile
    File sumstats
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
    Int window
    Int max_region_width
    Float window_shrink_ratio
    # can be helpful if adding finemapping with relaxed threshold after more stringent has already ben run.
    # does not include regions with lead snp < this
    Float p_threshold
    Float? minimum_pval
    String? set_variant_id_map_chr

    command <<<

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
            for(i = 1; i <= NF; i++) {
                h[$i] = i
            }
            exists=ph in h
            if (!exists) {
                print "Phenotype:"ph" not found in the given phenotype file." > "/dev/stderr"
                err = 1
                exit 1
            }
        }
        NR > 1 && $h[ph] != "NA" {
            vals[$h[ph] + 0] += 1
            print $1 > ph".incl"
            if ($h[ph] != 0 && $h[ph] != 1 && !err) {
                print "Phenotype:"ph" seems a quantitative trait. Setting var_y = 1 and prior_std = 0.05." > "/dev/stderr"
                print 1.0 > "var_y.txt"
                print 0.05 > "prior_std.txt"
                err = 1
            }
        }
        END {
            if (!err) {
                phi = vals["1"] / (vals["1"]+vals["0"])
                var_y = phi * (1-phi)
                std = 0.05 * sqrt(phi*(1-phi))
                print var_y > "var_y.txt"
                print std > "prior_std.txt"
            }
        }'

        if [[ $? -ne 0 ]]
        then
            echo "Error occurred while getting case control counts for ${pheno}"
            exit 1
        fi

        wc -l ${pheno}.incl | cut -f1 -d' ' > n_samples.txt

        make_finemap_inputs.py \
            --sumstats ${sumstats} \
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
            --exclude-MHC \
            --no-upload \
            --prefix ${pheno} \
            --out ${pheno} \
            --window ${window} \
            --max-region-width ${max_region_width} \
            --window-shrink-ratio ${window_shrink_ratio} \
            ${true='--scale-se-by-pval ' false=' ' scale_se_by_pval} \
            ${true='--x-chromosome' false=' ' x_chromosome} \
            ${true='--set-variant-id ' false=' ' set_variant_id} \
            ${true='--set-variant-id-map-chr ' false=' ' defined(set_variant_id_map_chr)}${set_variant_id_map_chr} \
            --p-threshold ${p_threshold} \
            ${true='--min-p-threshold ' false='' defined(minimum_pval)}${minimum_pval} \
            --wdl

            res=`cat ${pheno}_had_results`

            if [ "$res" == "False" ]; then
                touch ${pheno}".z"
                touch ${pheno}".lead_snps.txt"
                touch ${pheno}".bed"
            fi
    >>>

    output {

        Int n_samples = read_int("n_samples.txt")
        Float prior_std = read_float("prior_std.txt")
        Float var_y = read_float("var_y.txt")
        File incl = pheno + ".incl"
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

task filter {

    File variant_file
    File sumstat
    String base = basename(sumstat,".gz")

    command <<<

        python3 - <<EOF > ${base}

        import sys
        import gzip
        variant_file="${variant_file}"
        sumstat="${sumstat}"
        variants = {}
        with gzip.open(variant_file, 'rt') as f:
            for line in f:
                variants[line.strip().replace('chr', '').replace('X', '23').replace('Y', '24').replace('MT', '25').replace('M', '25')] = True
        print("${variant_file}",file=sys.stderr)
        print("${sumstat}",file=sys.stderr)
        print("${base}",file=sys.stderr)
        with gzip.open(sumstat, 'rt') as f:
            l = f.readline().strip()
            print(l)
            print(l,file=sys.stderr)
            for line in f:
                line = line.strip()
                s = line.split('\t')
                chr = s[0].replace('chr', '').replace('X', '23').replace('Y', '24').replace('MT', '25').replace('M', '25')
                id = chr + ':' + s[1] + ':' + s[2] + ':' + s[3]
                if id in variants:
                    print(line)
        EOF
        bgzip ${base}
        tabix -s1 -b2 -e2 ${base}.gz

    >>>

    output {
        File out = base+".gz"
        File out_tbi = base + ".gz.tbi"
    }

    runtime {

        docker: "eu.gcr.io/finngen-refinery-dev/bioinformatics:0.7"
        cpu: 1
        # 40M variants in variant_file to look up takes about 4G
        memory: "4 GB"
        disks: "local-disk 100 HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}

workflow finemap {

    String zones
    String docker
    String sumstats_pattern
    File phenolistfile
    File phenotypes

    Array[String] phenos = read_lines(phenolistfile)
    String? set_variant_id_map_chr

    scatter (pheno in phenos) {

        File sumstats = sub(sumstats_pattern,"\\{PHENO\\}",pheno)
        call filter{
            input: sumstat = sumstats
        }

        call preprocess {
            input: zones=zones, docker=docker, pheno=pheno, phenofile=phenotypes,
                sumstats=filter.out,set_variant_id_map_chr=set_variant_id_map_chr
        }

        if(preprocess.had_results) {
            call sub.ldstore_finemap {
                input: zones=zones, docker=docker, pheno=pheno,
                    n_samples=preprocess.n_samples, prior_std=preprocess.prior_std, var_y=preprocess.var_y,
                    incl=preprocess.incl, zfiles=preprocess.zfiles,
                    pheno=pheno, set_variant_id_map_chr=set_variant_id_map_chr
            }
        }
    }

    output {
        Array[File] bed = preprocess.bed 
        Array[Boolean] had_results = preprocess.had_results
        Array[File] out_susie_snp_filtered = select_all(ldstore_finemap.out_susie_snp_filtered)
        Array[File] out_susie_cred_summary = select_all(ldstore_finemap.out_susie_cred_summary)
        Array[File] out_susie_snp_filtered_99 = select_all(ldstore_finemap.out_susie_snp_filtered_99)
        Array[File] out_susie_cred_summary_99 = select_all(ldstore_finemap.out_susie_cred_summary_99)
        Array[File] out_susie_snp_filtered_extend = select_all(ldstore_finemap.out_susie_snp_filtered_extend)
        Array[File] out_susie_cred_summary_extend = select_all(ldstore_finemap.out_susie_cred_summary_extend)
        Array[File] out_susie_snp = select_all(ldstore_finemap.out_susie_snp)
        Array[File] out_susie_snp_tbi = select_all(ldstore_finemap.out_susie_snp_tbi)
        Array[File] out_susie_cred = select_all(ldstore_finemap.out_susie_cred)
        Array[File] out_susie_cred_99 = select_all(ldstore_finemap.out_susie_cred_99)

        Array[Array[File]] out_susie_rds = select_all(ldstore_finemap.out_susie_rds)
        Array[Array[Array[File]]] out_finemap_cred_regions = select_all(ldstore_finemap.finemap_cred_regions)
        Array[File] out_finemap_snp = select_all(ldstore_finemap.out_finemap_snp)
        Array[File] out_finemap_snp_tbi = select_all(ldstore_finemap.out_finemap_snp_tbi)
        Array[File] out_finemap_config = select_all(ldstore_finemap.out_finemap_config)
        Array[File] out_finemap_region = select_all(ldstore_finemap.out_finemap_region)
    }
}
