#!/usr/bin/env python2
import argparse
from logging import getLogger, StreamHandler, FileHandler, Formatter, DEBUG
from subprocess import Popen,PIPE, check_output
import pandas as pd
import gzip
from collections import defaultdict, namedtuple
import operator
import re
import tempfile
logger = getLogger(__name__)

log_format = Formatter('%(message)s')
handler = StreamHandler()
handler.setFormatter(log_format)
handler.setLevel(DEBUG)
logger.setLevel(DEBUG)
logger.addHandler(handler)
logger.propagate = False

class InvalidInputException(Exception):
    pass


class Region:
    def __init__(self,chr,start,stop):
        self.chr = chr.lstrip("chr")
        self.start=start
        self.stop=stop

class CSReport:
    def __init__(self,cs_id,comp_col,comp_val, comp_func, best_row, cs_dat):
        self.cs_id = cs_id
        self.comp_col = comp_col
        self.comp_val = comp_val
        self.comp_func = comp_func
        self.best_row = best_row
        self.cs_dat = cs_dat

def process_cred(cred_file, min_r2):
    REQ_COLUMNS=["trait","region","cs","cs_log10bf","cs_avg_r2","cs_min_r2","cs_size"]
    cred = pd.read_csv(cred_file,compression="gzip", sep="\t")
    if not all([col in cred.columns for col in REQ_COLUMNS]):
        raise InvalidInputException("All required columns not present in credible set file.")

    cred["good_cs"] = cred.cs_min_r2 > min_r2
    cred["cs_id"] = cred.apply( lambda row: row.region + str(row.cs),
                                axis=1)
    outcols = cred.columns.tolist()
    return (outcols,{t.cs_id:t for t in cred.itertuples(index=False)})


def get_variant_annots( regions, annot_file, outcols=["gene_most_severe","most_severe"], cpra=["chr","pos","ref","alt"] ):
    regions_file = tempfile.NamedTemporaryFile('w')
    for r in regions:
        regions_file.write('{}\t{}\t{}'.format(r.chr, r.start, r.stop) + "\n")

    header = ["chr","pos","ref","alt","gene_most_severe","most_severe"]
    hi = { h:i for i,h in enumerate(header)}

    miscols = [ c for c in cpra + outcols if c not in hi ]
    if len(miscols)>0:
        raise InvalidInputException("All required columns not in annotation file. Missing:" + ",".join(miscols))
    va = {}
    regions_file.flush()
    p = Popen(["tabix","-R",regions_file.name,annot_file ],stdout=PIPE, stderr=PIPE)

    for v in p.stdout:
        vd = v.strip().split("\t")
        varid="{}:{}:{}:{}".format( vd[hi[cpra[0]]], vd[hi[cpra[1]]], vd[hi[cpra[2]]], vd[hi[cpra[3]]])
        va[varid] = [ vd[hi[c]] for c in outcols ]

    return va

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('susie_cred', type=str)
    parser.add_argument('susie_snp', type=str)
    parser.add_argument('outprefix', type=str)

    parser.add_argument('--min_r2', type=float, default=0.25)
    parser.add_argument('--snp_outcols', type=str,
                    default="trait,region,v,cs,cs_specific_prob,chromosome,position,allele1,allele2,maf,beta,p,se")

    parser.add_argument('--cs_sum_snp_cols', type=str,
                    default="v,rsid,p,beta,sd,prob,cs,cs_specific_prob")

    parser.add_argument('--variant_annot', type=str)
    parser.add_argument('--variant_annot_cols', type=str)

    args = parser.parse_args()
    creds = process_cred(args.susie_cred, args.min_r2)
    best_vars = {}

    regions_of_cs=[]
    split_re = re.compile('[:-]')
    for cs_id, csdef in creds[1].iteritems():
        comp_col = "cs_specific_prob" if csdef.good_cs else "p"
        def smaller(a,b):
             return a is None or operator.lt(b,a)
        def bigger(a,b):
            return a is None or operator.gt(b,a)
        comp_func =  bigger if csdef.good_cs else smaller
        best_vars[cs_id] = CSReport(cs_id=cs_id, comp_col=comp_col, comp_val=None,
                                            comp_func=comp_func, best_row=[], cs_dat=csdef)

        regions_of_cs.append(Region(*split_re.split(csdef.region)))

    sum_snp_cols = [ c.strip() for c in args.cs_sum_snp_cols.split(",")]
    snp_cols = [c.strip() for c in args.snp_outcols.split(",")]

    var_annot_cols=[]
    var_annot={}
    if args.variant_annot and args.variant_annot_cols:
        var_annot_cols = args.variant_annot_cols.strip().split(",")
        var_annot = get_variant_annots(regions_of_cs,args.variant_annot, var_annot_cols)

    with open(args.outprefix + ".snp.filter.tsv", 'wt') as snpfile:
        with gzip.open(args.susie_snp, mode='rt') as snps:
            hl = snps.readline()
            hldat = hl.rstrip("\n").split("\t")
            hi = {h:i for i,h in enumerate(hldat)}
            hlen = len(hi)

            missing_sum_snp = set([ c for c in sum_snp_cols + snp_cols if c not in hi ])

            if len(missing_sum_snp)>0:
                raise InvalidInputException("All specified columns (missing:" + ",".join(missing_sum_snp)
                    + ") were not in snp file")

            snpfile.write("\t".join([hldat[hi[col]] for col in snp_cols] + var_annot_cols)+ "\n")
            for sl in snps:
                if sl == "":
                    continue
                ldat = sl.rstrip("\n").split("\t")
                if ldat[hi["cs"]] == "-1":
                    continue

                if(len(ldat) != hlen):
                    logger.warning("Wrong number of columns in row:" + sl)

                cs_id = ldat[hi["region"]] + ldat[hi["cs"]]
                cred_defs = best_vars[cs_id]

                varid = ldat[hi["v"]]
                va = [] if len(var_annot_cols)==0 else var_annot.get( varid, ["NA"] * len(var_annot_cols) )

                if cred_defs.cs_dat.good_cs:
                    snpfile.write("\t".join( [ldat[hi[col]] for col in snp_cols] + va) + "\n")

                compval = float(ldat[hi[cred_defs.comp_col]])
                if cred_defs.comp_func( cred_defs.comp_val, compval ):
                    cred_defs.best_row = ldat
                    cred_defs.comp_val = float(compval)

    with open(args.outprefix + ".cred.summary.tsv", mode='wt') as sumfile:
        sumfile.write( "\t".join(creds[0]) + "\t" + "\t".join([ hldat[hi[col]] for col in sum_snp_cols ] + var_annot_cols ) + "\n" )
        for cs_id, rep in best_vars.iteritems():
            orig_cs_dat = creds[1][cs_id]
            varid = rep.best_row[hi["v"]]
            va = var_annot.get( varid, ["NA"] * len(var_annot_cols))
            sumfile.write("\t".join(  [ str(d) for d in list(orig_cs_dat) ] )  + "\t" +
                            "\t".join([  rep.best_row[hi[col]] for col in sum_snp_cols ] + va ) + "\n")
