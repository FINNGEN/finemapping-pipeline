#!/usr/bin/env python3
import shlex
import sys
import requests
from typing import Dict
import argparse
import functools
import pysam
import tempfile
import subprocess
import io
import time
from math import ceil

import multiprocessing
from typing import List,Tuple

from collections import namedtuple

def get_ld_vars( chrom:str, pos:int, ref:str, alt:str, r2:float, ld_w:int, ld_source:str, retries:int=5):
    snooze=2
    snooze_multiplier = 2
    max_snooze=256

    ld_w=max(min(5000000,ld_w), 500000)
    url = "http://api.finngen.fi/api/ld?variant={chrom}:{pos}:{ref}:{alt}&panel={ld_source}&window={ld_w}&r2_thresh={r2}".format(chrom=chrom,
    pos=pos, ref=ref, alt=alt, ld_source=ld_source, ld_w=ld_w, r2=r2)

    print("Requestion {}".format(url) ,file=sys.stderr)
    sys.stderr.flush()
    r = requests.get(url)
    retries_left = retries
    while r.status_code!=200 and retries_left>0:
        retries_left-=1
        if r.status_code!=200:
            print("Error requesting ld for url {}. Error code: {}".format(url, r.status_code) ,file=sys.stderr)
            time.sleep(snooze)
            snooze = min(snooze * snooze_multiplier, max_snooze)
        r = requests.get(url)

    if r.status_code!=200:
        raise Exception("LD server response failed for {} attempts".format(retries) )

    return(r.json()["ld"])

def get_ld_vars_tomahawk(chrom:str, pos:int, ref:str, alt:str, r2:float, ld_w:int,
tomahawk_file_pattern_map:Dict[str,str], tomahaw_map_tabix, n_cpu=1):
    cpra = "{}:{}:{}:{}".format(chrom, pos, ref, alt)
    mapping = _get_tomahawk_region_mapping(cpra,tomahaw_map_tabix, ld_w)

    print("Getting ldvars for {}".format(cpra),file=sys.stderr)
    sys.stderr.flush()
    twkpos = mapping[1]

    twkcp = mapping[0].split(':')
    filename = tomahawk_file_pattern_map[chrom]
    tmp = next(tempfile._get_candidate_names())
    try:
        cmd_scalc = "tomahawk scalc -i {} -o {} -I {}:{}-{} -w {} -t {}".format(
        filename,tmp, twkcp[0], int(twkcp[1])-1, twkcp[1], ld_w, n_cpu)
        ret = subprocess.run(shlex.split(cmd_scalc), stdout=subprocess.DEVNULL)
        ret.check_returncode()

        view_scalc = "tomahawk view -i {}.two -H -r {}".format(tmp,r2)
        sb = subprocess.run(shlex.split(view_scalc), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,  check=True)

        ldres = parse_ld(io.TextIOWrapper(io.BytesIO(sb.stdout), encoding=sys.stdout.encoding), cpra, twkpos)

    finally:
       subprocess.call(['rm', tmp + ".two"])
    return ldres

def parse_ld(data_io, cpra, twk2cpra):
    """
    Parses the given tomahawk LD output using the given query variant and tomahawk_position-to-variant mapping
    Returns a list of dicts with keys variation1,variation2,r2,d_prime where variation1 is the query variant
    """
    for line in data_io:
        if not line.startswith('#') and line != '':
            break
    hdr = {h:i for i,h in enumerate(line.strip().split('\t'))}
    res = []
    used = {}
    for line in data_io:
        s = line.strip().split('\t')
        var1 = s[hdr['CHROM_A']] + ':' + s[hdr['POS_A']]
        var2 = s[hdr['CHROM_B']] + ':' + s[hdr['POS_B']]
        if var1 not in twk2cpra:
            app.logger.warning(var1 + ' tomahawk position not in given mapping, this should not happen. Ignoring')
        elif var2 not in twk2cpra:
            app.logger.warning(var2 + ' tomahawk position not in given mapping, this should not happen. Ignoring')
        else:
            var1 = twk2cpra[var1]
            var2 = twk2cpra[var2]
            if var2 == cpra:
                temp = var1
                var1 = var2
                var2 = temp
            if var1 == cpra and var2 not in used:
                res.append({'variation1': var1, 'variation2': var2, 'r2': round(float(s[hdr['R2']]), 2), 'd_prime': round(float(s[hdr['DPrime']]), 2) } )
                used[var2] = True
    return res

def _get_tomahawk_region_mapping(cpra, tomahaw_map_tabix, window):
    """
    Gets tomahawk position for the query variant and a tomahawk_position-to-variant mapping for variants in the given panel within the given window
    Returns a tuple (tomahawk chr:pos, dict from tomahawk chr:pos to chr:pos:ref:alt)
    Raises if variant not found
    """
    twk2cpra = {} # mapping from tomahawk positions chr:pos to actual variants chr:pos:ref:alt
    twk = None # tomahawk position of query variant
    s = cpra.split(':')

    tabix_iter = tomahaw_map_tabix.fetch(s[0], max(1,int(s[1])-round(window)-1), int(s[1])+round(window), parser=None)

    for row in tabix_iter:
        s = row.split('\t')
        if s[2] == cpra:
            twk = s[3]
        twk2cpra[s[3]] = s[2]
    if twk is None:
        raise Exception('Variant {} not found: '.format(cpra))
    return (twk, twk2cpra)


def get_span( ldlist, from_bp ):

    min=int(from_bp)
    max=int(from_bp)
    for ld in ldlist:
        var = ld["variation2"].split(":")
        bp = int(var[1])
        if bp>max:
            max=int(var[1])
        if bp<min:
            min=int(var[1])

    return (min,max)


VariantData = namedtuple('VariantData',['chrom','pos','ref','alt','pval'])

def process_variants( varlist:List[VariantData], r2, ld_dao, ld_width:int, min_buffer, disable_pvalspeedup=False) -> List[Tuple[str,int,int]]:
    regions=[]

    if not disable_pvalspeedup:
        varlist = sorted(varlist, key=lambda x: x.pval)

    done={}

    for v in varlist:

        vid = ":".join([ v.chrom, str(v.pos), v.ref, v.alt])
        if (not disable_pvalspeedup) and vid in done:
            ## more significant hit as already LD eaten this, speedup by not needing to compute LD.
            print("skipping {} as already eaten by stronger variant {}".format(vid, done[vid]), file=sys.stderr )
            continue

        ldpairs = ld_dao(v.chrom,v.pos,v.ref,v.alt, r2=r2,
        ld_w=ld_width)

        if not disable_pvalspeedup:
            for p in ldpairs:
                done[ p["variation2"] ] = p["variation1"]

        span=get_span(ldpairs, from_bp=v.pos)

        left = min(span[0], v.pos- min_buffer)
        right = max(span[1], v.pos + min_buffer)

        regions.append([v.chrom,left, right])

    return _merge_regions(regions)

def _merge_regions(regions):
    regions = sorted(regions, key=lambda x: (x[0],x[1]))
    merged = [regions[0]]
    for i in range(1,len(regions)):
        curr = regions[i]
        if curr[0]==merged[-1][0] and curr[1]<=merged[-1][2]:
            merged[-1][2]=max(curr[2], merged[-1][2])
        else:
            merged.append( curr )
    return merged


def _backforth(func, queue, *args):
    res = func(*args)
    queue.put(res)

def main(args):
    regions=list()

    varlist= []

    ## if pvalue given in second column, does speedup by not getting ld for variants already covered by stronger variants
    disable_pvalspeedup = False

    for v in sys.stdin:
        v = v.strip().split("\t")
        cpra = v[0].split("_")
        pval= None
        if len(v)>1:
            pval=float(v[1])
        else:
            disable_pvalspeedup=True

        varlist.append(VariantData(cpra[0],int(cpra[1]),cpra[2],cpra[3],pval))


    varlist = sorted(varlist, key=lambda x:(x[0],x[1]))

    chunksize =  ceil(len(varlist)/args.n_cpu)
    chunks =[ varlist[i*(chunksize):min(len(varlist),(i*chunksize+chunksize))] for i in range(0, args.n_cpu) ]
    queue = multiprocessing.Queue()
    processes = []
    qout = multiprocessing.Queue()

    print("pvalue speedup disabled {}".format(disable_pvalspeedup), file=sys.stderr)
    for c in chunks:
        if  args.ld_source is None:
            filepattern_map = { str(chr):args.tomahawk_pattern.replace("[CHROMOSOME]",str(chr)) for chr in range(1,24) }
            ###  TODO check  if all files exist in initialization
            ## ld source
            tomahaw_map_tabix =  pysam.TabixFile(args.tomahawk_mapping, parser=None)
            ld_fetch = functools.partial( get_ld_vars_tomahawk, tomahawk_file_pattern_map=filepattern_map,
            tomahaw_map_tabix=tomahaw_map_tabix, n_cpu=args.tomahawk_cpu)
        else:
            ld_fetch = functools.partial(get_ld_vars, ld_source=args.ld_source)
        processes.append( multiprocessing.Process(target=_backforth, args=(process_variants, qout, c, args.r2, ld_fetch, args.search_width, args.min_buffer, disable_pvalspeedup)  ))

    for p in processes:
        p.start()

    finished=0
    while finished<len(processes):
        for i in range(len(processes)-1,-1,-1):
            if processes[i].exitcode and  processes[i].exitcode!=0:
                [ p.terminate() for p in processes]
                raise Exception("Exception occurred in subprocess. Killing all.")
            elif processes[i].exitcode==0:
                finished+=1
        time.sleep(1)

    results=[]
    for i in range(len(processes)):
        result = qout.get()
        results.extend(result)

    final = _merge_regions(results)
    print("\n".join([ "\t".join([ str(v) for v in  f])  for f in final ]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--r2', default=0.6, help="r2 threshold")
    parser.add_argument('--search_width', type=int,default=5000000, help="ld extent to search LD partners")
    parser.add_argument('--min_buffer', default=20000, help="Buffer basepairs from variant")
    parser.add_argument('--ld_source', help="LD panel if using LD server (sisu3,sisu4). If given then LD server will be used")
    parser.add_argument('--n_cpu', default=1, type=int,  help="Splits the input variants to n_cpu chunks and runs in parallel.")
    parser.add_argument('--tomahawk_pattern', default="/mnt/ld/sisu4/chr[CHROMOSOME]_phased_SNPID.twk", help="LD panel")
    parser.add_argument('--tomahawk_mapping', default="/mnt/ld/sisu4/sisu4_twk_mapping.txt.gz", help="LD panel")
    parser.add_argument('--tomahawk_cpu', type=int, default=1)

    args = parser.parse_args()
    main(args)
