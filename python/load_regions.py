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