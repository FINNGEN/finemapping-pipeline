path=$1
target=$2
suffixes=("**/*.config.bgz" "**/*region.bgz" "**/*snp.bgz" "**/*cred.bgz" "**/*.cred*" "**/*.config" "**/*.snp")


for suf in ${suffixes[*]} 
do
    echo "gsutil -m cp"$path"/"$suf $target
done


