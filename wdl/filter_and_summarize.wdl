
task do_filt_summ {
    String zone
    String docker
    Int cpu
    Int mem
    String pheno
    File susie_snps
    File susie_snps_tbi
    File susie_cred
    File susie_cred_99
    File? snp_annot_file
    File? snp_annot_file_tbi
    String? snp_annotation_fields
    String? set_variant_id_map_chr
    Float good_cred_r2

    command <<<

        filter_and_summarize.py --min_r2 0 ${susie_cred} \
            ${susie_snps} ${pheno}.SUSIE_all_cred \
            ${true='--variant_annot ' false='' defined(snp_annot_file)}${snp_annot_file} \
            ${true='--variant_annot_cols ' false='' defined(snp_annotation_fields)}${snp_annotation_fields} \
            ${true='--set_variant_id_map_chr ' false='' defined(set_variant_id_map_chr)}${set_variant_id_map_chr}

        filter_and_summarize.py --min_r2 ${good_cred_r2} ${susie_cred} \
            ${susie_snps} ${pheno}.SUSIE \
            ${true='--variant_annot ' false='' defined(snp_annot_file)}${snp_annot_file} \
            ${true='--variant_annot_cols ' false='' defined(snp_annotation_fields)}${snp_annotation_fields} \
            ${true='--set_variant_id_map_chr ' false='' defined(set_variant_id_map_chr)}${set_variant_id_map_chr}

        filter_and_summarize.py --min_r2 ${good_cred_r2} ${susie_cred_99} \
            ${susie_snps} ${pheno}.SUSIE_99 \
            --cs_column "cs_99" --pip_column "cs_specific_prob_99" \
            --snp_outcols "trait,region,v,cs_99,cs_specific_prob_99,chromosome,position,allele1,allele2,maf,beta,p,se" \
            --snp_outcols_header "trait,region,v,cs,cs_specific_prob,chromosome,position,allele1,allele2,maf,beta,p,se" \
            --cs_sum_snp_cols "v,rsid,p,beta,sd,prob,cs_99,cs_specific_prob_99" \
            --cs_sum_snp_cols_header "v,rsid,p,beta,sd,prob,cs,cs_specific_prob" \
            ${true='--variant_annot ' false='' defined(snp_annot_file)}${snp_annot_file} \
            ${true='--variant_annot_cols ' false='' defined(snp_annotation_fields)}${snp_annotation_fields} \
            ${true='--set_variant_id_map_chr ' false='' defined(set_variant_id_map_chr)}${set_variant_id_map_chr}

        ## add 95% cs + extended snps from 99.
        filter_and_summarize.py --min_r2 ${good_cred_r2} ${susie_cred} \
            ${susie_snps} ${pheno}.SUSIE_extend \
            --extend_cs_column "cs_99" --extend_pip_column "cs_specific_prob_99" \
            --snp_outcols "trait,region,v,cs,cs_specific_prob,cs_99,cs_specific_prob_99,chromosome,position,allele1,allele2,maf,beta,p,se" \
            ${true='--variant_annot ' false='' defined(snp_annot_file)}${snp_annot_file} \
            ${true='--variant_annot_cols ' false='' defined(snp_annotation_fields)}${snp_annotation_fields} \
            ${true='--set_variant_id_map_chr ' false='' defined(set_variant_id_map_chr)}${set_variant_id_map_chr}

    >>>

    output {
        File out_susie_snp_filtered = pheno + ".SUSIE.snp.filter.tsv"
        File out_susie_cred_summary = pheno + ".SUSIE.cred.summary.tsv"
        File out_susie_snp_filtered_99 = pheno + ".SUSIE_99.snp.filter.tsv"
        File out_susie_cred_summary_99 = pheno + ".SUSIE_99.cred.summary.tsv"
        File out_susie_snp_filtered_extend = pheno + ".SUSIE_extend.snp.filter.tsv"
        File out_susie_cred_summary_extend = pheno + ".SUSIE_extend.cred.summary.tsv"
    }

    runtime {
        docker: "${docker}"
        cpu: "${cpu}"
        memory: "${mem} GB"
        disks: "local-disk 10 HDD"
        zones: "${zone}"
        preemptible: 2
        noAddress: true
    }
}

workflow filter_and_summarize {
    String res_snpfile_list
    Array[String] res_snpfiles = read_lines(res_snpfile_list)
    String docker
    String? set_variant_id_map_chr
    String zone
    Float good_cred_r2

    scatter (res in res_snpfiles) {
        String stem=sub(res,".SUSIE.snp.bgz$","")
        String pheno=basename(stem)
        call do_filt_summ {
            input: zone=zone, pheno=pheno, docker=docker, susie_snps=stem+".SUSIE.snp.bgz",
                susie_snps_tbi=stem+".SUSIE.snp.bgz.tbi", susie_cred=stem+".SUSIE.cred.bgz",
                susie_cred_99=stem+".SUSIE.cred_99.bgz", set_variant_id_map_chr=set_variant_id_map_chr,
                good_cred_r2=good_cred_r2
        }
    }
}
