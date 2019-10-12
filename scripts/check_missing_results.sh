#!/usr/bin/env bash

##
# utility script to check if all results are successfully copied to target bucket.
# NOTE THAT cred1 file can be missing where it's not output by finemap (not likely conf).
##
cromwell_preprocess_path=$1 ## e.g gs://fg-cromwell/finemap/6c939688-aff7-4340-a8f4-f8eaf2c1b3a2/call-preprocess/**/*.z
target_bucket=$2 #"gs://r4_data_west1/finemap/*"
out=$3
#gs://r4_data_west1/finemap/C3_DIGESTIVE_ORGANS.SUSIE.snp.bgz
#gs://r4_data_west1/finemap/C3_DIGESTIVE_ORGANS.chr13.108495450-111495450.config
#gs://r4_data_west1/finemap/C3_DIGESTIVE_ORGANS.chr13.108495450-111495450.cred1
#gs://r4_data_west1/finemap/C3_DIGESTIVE_ORGANS.chr13.108495450-111495450.cred2
#gs://r4_data_west1/finemap/C3_DIGESTIVE_ORGANS.chr13.108495450-111495450.snp
#gs://r4_data_west1/finemap/C3_DIGESTIVE_ORGANS.chr13.108495450-111495450.susie.cred
#gs://r4_data_west1/finemap/C3_DIGESTIVE_ORGANS.chr13.108495450-111495450.susie.snp

suffixes=(.cred1 .config .snp .susie.cred .susie.snp)
join -v 1 -1 1 -2 1 <(gsutil ls $cromwell_preprocess_path  | while read f; do n=`basename $f`;if [[ $n =~ .*[0-9]+-[0-9]+\.z$ ]]; then for suf in "${suffixes[@]}"; do  echo ${n/.z/}$suf; done; fi; done | sort -b -k 1,1)  <(gsutil ls $target_bucket  | while read f; do n=`basename $f`; echo $n; done | sort -b -k 1,1) >$out
