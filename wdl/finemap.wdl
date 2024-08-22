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
    File? manual_regions

    command <<<
cat << "__EOF__" > manual_regions.py
"""
Create zfiles with region definitions
"""
#import pandas as pd
#import numpy as np
import argparse,gzip,math
from typing import NamedTuple, List, Dict,IO,Iterator,Callable,Union,Optional
from collections import defaultdict
from contextlib import contextmanager
ZFILE_COLUMNS=['rsid', 'chromosome', 'position', 'allele1', 'allele2', 'maf', 'beta', 'se', 'p']

Datatype = Union[int,str,float,None]

class Region(NamedTuple):
    chrom: str
    start: int
    end: int

    def overlaps(self, other: 'Region')->bool:
        """Check if two Regions overlap
        """
        #check that both regions are valid
        if self.end < self.start:
            raise Exception(f"Region {self} is invalid: start is larger than end!")
        if other.end < other.start:
            raise Exception(f"Region {other} is invalid: start is larger than end!")
        if self.chrom == other.chrom:
            if (self.start <= other.end) and (other.start <= self.end):
                return True
        return False


def load_regions(regionfile: str, exclude_region: Region) -> List[Region]:
    """
    Load regions from a regionfile. regions are in form chr:pos-pos
    """
    out=[]
    with open(regionfile,"r") as f:
        for line in f:
            d = line.split(":")
            chrom = d[0]
            start = d[1].split("-")[0]
            end = d[1].split("-")[1]
            r = Region( str(chrom),int(start),int(end) )
            if not r.overlaps(exclude_region):
                out.append(r)
    return out

@contextmanager
def uopen(fname:str,oper_type:str)->Iterator[IO]:
    """
    Universal opener
    Open both gzipped and plaintext files with no fuzz
    """
    gz_magicnumber=b"\x1f\x8b"
    type="normal"
    with open(fname,"rb") as f:
        if f.read()[0:2] == gz_magicnumber:
            type="gz"
    if type=="normal":
        with open(fname,oper_type,encoding="utf-8") as f:
            yield f
    elif type == "gz":
        with gzip.open(fname,oper_type,encoding="utf-8") as f:
            yield f
    else:
        raise Exception("invalid file format")


def load_data(fname:str,columns:List[str],dtypes:Dict[str,Callable[[str],Datatype]])->Dict[str,List[Dict[str,Datatype] ]]:
    """
    Load data from possibly gzip-compressed file.
    Returns: Dictionary with chromosome as key, list of variants in that chromosome as the value.
    """
    output = defaultdict(list)
    with uopen(fname,"rt") as f:
        #read header
        header =f.readline().strip().split("\t")
        hdi = {a:i for i,a in enumerate(header)}
        try:
            _ = [hdi[a] for a in columns]
        except:
            raise Exception(f"Header did not contain columns that were defined: {[a for a in columns if a not in header]}")
        
        for l in f:
            ld = l.strip().split("\t")
            line_data = {}
            for c in columns:
                line_data[c] = dtypes[c](ld[hdi[c]])
            output[ld[hdi[columns[0]]]].append(line_data)
    return output

def tryfloat(value:str)->Optional[float]:
    try:
        return float(value)
    except:
        return None

def data_formatter(value:Datatype)->str:
    if value is None:
        return "NA"
    elif type(value) == int:
        return str(value)
    elif type(value) == float:
        if math.isnan(value):
            return "NA"
        return f"{value:.6g}"
    elif type(value) == str:
        return value
    else:
        raise Exception(f"unsupported type of value: {type(value)} for value {value}")

def extract_regions(summstat_file:str, pheno:str, regions:List[Region], cs:List[str], replacevals:Dict[str,str]):
    """
    Extract regions with no filtering. Write zfiles as soon as regions have been created.
    """
    if not regions:
        return
    coltypes = {
        cs[0]:str,   #chrom
        cs[1]:int,   #pos
        cs[2]:str,   #ref
        cs[3]:str,   #alt
        cs[4]:float, #maf
        cs[5]:tryfloat, #beta
        cs[6]:tryfloat, #se
        cs[7]:tryfloat  #pval
    }
    rename_d = {a:ZFILE_COLUMNS[i+1] for i,a in enumerate(cs)}
    print("starting data loading")
    data = load_data(summstat_file,cs,coltypes)

    print("starting region creation...")
    for r in regions:
        out_region = []
        for l in data[r.chrom]:
            if (l[cs[0]] == r.chrom) and (l[cs[1]] >= r.start) and (l[cs[1]] <= r.end):
                out_l = {}
                # chromosome
                # replace chrom with chrommap, add chr prefix
                out_l[rename_d[cs[0]]] = "chr"+ replacevals.get(l[cs[0]],l[cs[0]])
                # rsid
                out_l["rsid"] = f"{out_l[rename_d[cs[0]]]}_{l[cs[1]]}_{l[cs[2]]}_{l[cs[3]]}"
                # maf
                out_l[rename_d[cs[4]]] = min(l[cs[4]],1.0-l[cs[4]])

                #rest
                for c in cs:
                    if c not in [cs[0],cs[4]]:
                        out_l[rename_d[c]] = l[c]

                out_region.append(out_l)
        out_name = f"{pheno}.chr{r.chrom}.{r.start}-{r.end}.z"
        #write output
        with open(out_name,"wt") as f:
            header_line = " ".join(ZFILE_COLUMNS)+"\n"
            f.write(header_line)
            for ol in out_region:
                line = " ".join([data_formatter(ol[a]) for a in ZFILE_COLUMNS])+"\n"
                f.write(line)


if __name__ == "__main__":
    ap = argparse.ArgumentParser("Select regions from summstats, create zfiles")
    ap.add_argument("summstat")
    ap.add_argument("--phenotype",required=True)
    ap.add_argument("--region-file",required=True)
    ap.add_argument("--columns",required=True,nargs=8,help="columns (8 in total): chrom, pos, a1, a2, maf, beta, se, p")
    ap.add_argument("--rm-region",type=str,default=None,help="if a region overlaps this region, it is removed.")
    ap.add_argument("--chrmap",type=str)
    args=ap.parse_args()
    rm_region = Region("-1",-1,-1)
    if args.rm_region != None:
        rm_region = Region( args.rm_region.split(":")[0],int(args.rm_region.split(":")[1].split("-")[0]),int(args.rm_region.split(":")[1].split("-")[1]) )
        print(f"rm region: {rm_region}")
    regions = load_regions( args.region_file, rm_region)
    #replacevals = (args.chrmap.split("=")[0],args.chrmap.split("=")[1])
    replacevals = { m[0].strip():m[1].strip() for m in [x.split("=") for x in args.chrmap.split(",")] }
    print(f"regions: {regions}")
    extract_regions(args.summstat,args.phenotype, regions, args.columns,replacevals)
__EOF__
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

        #if custom region list, use that
        #else, normal region selection
        if ${true='true' false='false' defined(manual_regions)}; then
            python3 manual_regions.py ${sumstats} \
                --phenotype ${pheno} \
                --region-file ${manual_regions} \
                --columns "${chromosome_col}" ${position_col} ${allele1_col} ${allele2_col} ${freq_col} ${beta_col} ${se_col} ${p_col} \
                --rm-region 6:25000000-34000000\
                ${true='--chrmap ' false=' ' defined(set_variant_id_map_chr)} "${set_variant_id_map_chr}"
            echo "true" > ${pheno}_had_results
            touch ${pheno}".lead_snps.txt"
            touch ${pheno}".bed"
            echo "Manual regions created from supplied region file" >> tmp.log
            echo "The following regions were made:" >> tmp.log
            cat tmp.log ${manual_regions} > ${pheno}".log"
        else
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
