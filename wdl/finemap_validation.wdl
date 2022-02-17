version 1.0

workflow finemap_validation {
    input {
        Array[String] current
        Array[String] previous
        String docker
        String zones
    }
    call fmap_gather{
        input: cred_summary_file = current[1],release_tag = current[0],docker=docker,zones=zones
    }
    call fmap_gather as fmap_gather_previous{
        input: cred_summary_file = previous[1],release_tag = previous[0],docker=docker,zones=zones
    }
    call plot{
        input: current_files = [fmap_gather.pheno_level_summary,fmap_gather.cs_level_summary],
        previous_files = [fmap_gather_previous.pheno_level_summary,fmap_gather_previous.cs_level_summary],
        current_tag = current[0],previous_tag = previous[0],zones=zones ,docker=docker
    }

    output {
        Array[File] previous_files = [fmap_gather_previous.pheno_level_summary,fmap_gather_previous.cs_level_summary]
        Array[File] current_files = [fmap_gather.pheno_level_summary,fmap_gather.cs_level_summary]
        File plots = plot.plots
    }
}

task fmap_gather {
    input{
        String docker
        File cred_summary_file
        String release_tag
        String zones
    }
    Array[File] cred_summaries = read_lines(cred_summary_file)

    command <<<
        set -euxo pipefail
        cat > summarise.py <<EOF
        import os, sys
        #
        # Summarise data across all susie credset files
        # One: get cs/pheno. This one is calculated
        # Two: get distribution of variants in CS. This one can be just all of the data in the files, without headers. Easy to do even with awk or something.

        LIST_OF_FILES = sys.argv[1]
        CS_SUMMARY_NAME = sys.argv[2]+"_cs_per_pheno.tsv"
        cs_summ_header = ["phenotype","n_good_cs","n_bad_cs","n"]
        VAR_DIST_NAME = sys.argv[2]+"_cred_summaries.tsv"
        with open(LIST_OF_FILES) as f:
            files = [a.strip("\n") for a in f.readlines()]
        #get header
        with open(files[0]) as f:
            var_f_header = f.readline()
        var_file = open(VAR_DIST_NAME,"w")
        cs_file = open(CS_SUMMARY_NAME,"w")
        _=var_file.write(var_f_header)
        _=cs_file.write("\t".join(cs_summ_header)+"\n")
        
        for fname in files:
            with open(fname,"r") as infile:
                #var file takes all rows and takes the amount of variants in them, as well as maybe good_cs
                #cs file only counts the good and bad cs
                header = infile.readline().rstrip("\r\n").split("\t")
                h_idx = {a:i for i,a in enumerate(header)}
                n_good_cs = 0
                n_bad_cs = 1
                pheno=""
                for line in infile:
                    data = line.rstrip("\r\n").split("\t")
                    good_cs = data[h_idx["good_cs"]]
                    pheno = data[h_idx["trait"]]
                    if good_cs == "True":
                        n_good_cs += 1
                    else:
                        n_bad_cs += 1
                    #write cs data to file
                    _=var_file.write(line)
                cs_file_row = "\t".join( (pheno,str(n_good_cs),str(n_bad_cs),str(n_good_cs+n_bad_cs)) )+"\n"
                cs_file.write(cs_file_row)

        var_file.close()
        cs_file.close()
        EOF
        python3 summarise.py "~{write_lines(cred_summaries)}" "~{release_tag}"

    >>>

    output{
        File pheno_level_summary = release_tag + "_cs_per_pheno.tsv"
        File cs_level_summary = release_tag + "_cred_summaries.tsv"
    }

    runtime {
        docker: "${docker}"
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 20 HDD"
        zones: "${zones}"
        preemptible: 2
        noAddress: true
    }
}

task plot{
    input{
        Array[File] current_files
        Array[File] previous_files
        String current_tag
        String previous_tag
        String docker
        String zones

    }
    output {
        File plots = current_tag + "_" +previous_tag+ "_plots.pdf"
    }

    command <<<
        #script
        cat > script.R <<'EOF'
        library(data.table)
        library(ggplot2)
        library(tidyverse)

        previous_pheno <- fread("~{previous_files[0]}")
        previous_cs <- fread("~{previous_files[1]}")
        current_pheno <- fread("~{current_files[0]}")
        current_cs <- fread("~{current_files[1]}")
        current_tag <- "~{current_tag}"
        previous_tag <- "~{previous_tag}"
        #join the phenotype files
        pheno_join <- current_pheno %>% left_join(previous_pheno,by="phenotype",suffix = c(".curr",".prev"))


        #aggregate the sizes etc per pheno in the other file, join to ther
        #good cs only
        previous_agg_good <- previous_cs %>%filter(good_cs==TRUE)%>% group_by(trait) %>% summarise(avg_size = mean(cs_size))
        current_agg_good <- current_cs %>%filter(good_cs==TRUE) %>% group_by(trait) %>% summarise(avg_size = mean(cs_size))
        joined <- current_agg_good %>% left_join(previous_agg_good,by="trait",suffix = c(".curr",".prev"))
        #pdf ON
        pdf("~{current_tag}_~{previous_tag}_plots.pdf",onefile=TRUE)

        #n credsets per pheno prev vs current
        model <- lm(n.curr ~ 0 + n.prev,data=pheno_join)
        slope <- model$coeff[1]
        ggplot(pheno_join,aes(x=n.prev,y=n.curr))+geom_point()+
        geom_abline(intercept=0,slope=slope)+
        geom_abline(intercept=0,slope=1,linetype=5)+
        xlab(paste(previous_tag, "# of cs per phenotype"))+
        ylab(paste(current_tag, "# of cs per phenotype"))+
        ggtitle(paste("# of cs per phenotype,", previous_tag,"vs",current_tag))

        #lm
        model <- lm(avg_size.curr ~ 0 + avg_size.prev,data=joined)
        slope <- model$coeff[1]
        ggplot(joined,aes(x=avg_size.prev,y=avg_size.curr))+geom_point()+
        geom_abline(intercept=0,slope=slope)+
        geom_abline(intercept=0,slope=1,linetype=5)+
        xlab(paste(previous_tag,"average cs size"))+
        ylab(paste(current_tag,"average cs size"))+
        ggtitle(paste("Average cs size",previous_tag,"vs",current_tag))

        #plot cs size vs avg r2, min r2 to see what's what
        ggplot(previous_cs,aes(x=cs_avg_r2,y=cs_size,color=good_cs))+geom_point()+
        xlab(paste(previous_tag,"CS average r2"))+
        ylab(paste(previous_tag,"CS size"))+
        ggtitle(paste("Avg r2 vs size,",previous_tag))

        ggplot(current_cs,aes(x=cs_avg_r2,y=cs_size,color=good_cs))+geom_point()+
        xlab(paste(current_tag,"CS average r2"))+
        ylab(paste(current_tag,"CS size"))+
        ggtitle(paste("Avg r2 vs size,",current_tag))

        ggplot(previous_cs,aes(x=cs_min_r2,y=cs_size,color=good_cs))+geom_point()+
        xlab(paste(previous_tag,"CS min r2"))+
        ylab(paste(previous_tag,"CS size"))+
        ggtitle(paste("Min r2 vs size,",previous_tag))

        ggplot(current_cs,aes(x=cs_min_r2,y=cs_size,color=good_cs))+geom_point()+
        xlab(paste(current_tag,"CS min r2"))+
        ylab(paste(current_tag,"CS size"))+
        ggtitle(paste("Min r2 vs size,",current_tag))

        dev.off()
        EOF
        #for each of the inputs, run the script
        Rscript script.R
    >>>

    runtime {
        docker: "${docker}"
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 20 HDD"
        zones: "${zones}"
        preemptible: 2
        noAddress: true
    }
}

