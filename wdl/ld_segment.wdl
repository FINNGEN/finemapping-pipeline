

task ld_seq {

  File sumstat
  String chrcol = "#chrom"
  String poscol = "pos"
  String refcol = "ref"
  String altcol = "alt"
  String pcol = "pval"
  Array[File] tomahawks
  File thawk_map
  File thawk_map_tbi = thawk_map + ".tbi"
  File pycode = "gs://r9_data/ldseq/ld_segment.py"
  Float p_threshold = 0.00000005
  Float r2 = 0.4
  Int threads = 10
  Int tomahawk_cpu = 1

  String filepattern = "chr[CHROMOSOME]_phased_SNPID.twk"

  String bedfile = basename(sumstat) + ".bed"

  String docker= "eu.gcr.io/finngen-refinery-dev/tomahawk:beta-0.7.1-dirty-fg-v1"
  Int cpu = 8
  Int mem = 8
  Int disk=5000

  command <<<
      set -euxo pipefail

      echo ${sumstat}
      for f in ${sep=" " tomahawks}; do
        mv $f .
      done


      zcat ${sumstat} | awk 'BEGIN{ FS=OFS="\t"}
                            NR==1{ for(i=1;i<=NF;i++)  { h[$i]=i } }
                            NR>1 && $h["${pcol}"]<${p_threshold}{ print $h["${chrcol}"]"_"$h["${poscol}"]"_"$h["${refcol}"]"_"$h["${altcol}"],$h["${pcol}"] }
                            ' > hits

      head hits

      nhits=$(wc -l hits | awk '{ print $1}')

      chmod 777 ${pycode}

      echo "N hits $nhits"
      if [[ $nhits -lt 1 ]]; then
        echo "No significant variants identified."
        echo "False" > had_results
        touch ${bedfile}
      else
        cat hits | ld_segment.py --r2 ${r2} --tomahawk_mapping ${thawk_map} \
          --tomahawk_pattern ${filepattern}  --n_cpu 10 --tomahawk_cpu ${tomahawk_cpu} > ${bedfile}
          echo "True" > had_results
      fi
  >>>

  output {
      File segments = bedfile
      Boolean had_results = read_boolean("had_results")
  }


  runtime {
      docker: "${docker}"
      cpu: "${cpu}"
      memory: "${mem} GB"
      disks: "local-disk ${disk} HDD"
      preemptible: 2
      zones: "europe-west1-b europe-west1-c europe-west1-d"
      noAddress: true
  }


}



workflow ld_segment {

  String sumstatsfile
  Array[String] sumstats = read_lines(sumstatsfile)

  scatter (stat in sumstats) {
    call ld_seq { input: sumstat=stat }
  }

}
