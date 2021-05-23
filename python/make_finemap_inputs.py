#!/usr/bin/env python
import argparse
import os
import os.path
import json
import numpy as np
import scipy as sp
from scipy import stats
import pandas as pd
import pybedtools
from collections import OrderedDict, defaultdict
from logging import getLogger, StreamHandler, FileHandler, Formatter, DEBUG
from pybedtools import BedTool
from utils import run_command, run_command_bg, make_executable

# logger
logger = getLogger(__name__)
# log_format = Formatter('%(asctime)s - %(levelname)s - %(message)s')
log_format = Formatter('%(message)s')
handler = StreamHandler()
handler.setFormatter(log_format)
handler.setLevel(DEBUG)
logger.setLevel(DEBUG)
logger.addHandler(handler)
logger.propagate = False

CHROM_CONSTANT = int(1e11)
# mapping chromosome string (incluing 01-09 and X) and integer
CHROM_MAPPING_INT = dict([(str(i), i) for i in range(1, 24)] + [('0' + str(i), i) for i in range(1, 10)] + [('X', 23)])
# mapping chromosome string (incluing 01-09, 23, and X) and correct string
CHROM_MAPPING_STR = dict([(str(i), str(i)) for i in range(1, 23)] + [('0' + str(i), str(i)) for i in range(1, 10)] +
                         [('X', 'X'), ('23', 'X')])
FINEMAP_COLUMNS = ['rsid', 'chromosome', 'position', 'allele1', 'allele2', 'maf', 'beta', 'se', 'p']


def convert_chrpos(chromosome, position):
    return CHROM_CONSTANT * chromosome + position


def read_sumstats(path,
                  rsid_col='rsid',
                  chromosome_col='chromosome',
                  position_col='position',
                  allele1_col='allele1',
                  allele2_col='allele2',
                  maf_col='maf',
                  freq_col='freq',
                  beta_col='beta',
                  se_col='se',
                  flip_col='flip',
                  p_col='p',
                  delimiter='\s+',
                  recontig=False,
                  set_variant_id=False,
                  variant_id_chr_map={},
                  flip_beta=False,
                  grch38=False,
                  scale_se_by_pval=False,
                  extra_cols=None):
    logger.info("Loading sumstats: " + path)
    sumstats = pd.read_csv(path,
                           delimiter=delimiter,
                           dtype={
                               chromosome_col: str,
                               position_col: int
                           },
                           compression='gzip' if path.endswith('gz') else 'infer')

    req_cols = [chromosome_col, position_col, allele1_col, allele2_col, beta_col, se_col]
    if not set_variant_id:
        req_cols.append(rsid_col)

    if extra_cols is not None:
        req_cols = req_cols + extra_cols
    missing = [c for c in req_cols if c not in sumstats.columns]

    if len(missing) > 0:
        logger.error("All required columns not present in the data. Missing columns: " + " ".join(missing))
        raise Exception("All required columns not present in the data. Missing columns: " + " ".join(missing))

    sumstats = sumstats.rename(index=str,
                               columns={
                                   rsid_col: 'rsid',
                                   chromosome_col: 'chromosome',
                                   position_col: 'position',
                                   allele1_col: 'allele1',
                                   allele2_col: 'allele2',
                                   beta_col: 'beta',
                                   se_col: 'se'
                               })

    if grch38:
        sumstats['chromosome'] = sumstats.chromosome.str.replace('^chr', '')
    # remove non-compatible chromosomes
    sumstats = sumstats.loc[sumstats.chromosome.isin(CHROM_MAPPING_STR.keys()), :]
    # convert to chrpos: 1e11 * chr + pos
    chromosome_int = sumstats.chromosome.map(CHROM_MAPPING_INT).astype(int)
    sumstats['chrpos'] = convert_chrpos(chromosome_int, sumstats.position)

    if recontig:
        sumstats['chromosome'] = np.where(chromosome_int < 10, '0' + sumstats.chromosome, sumstats.chromosome)

    sumstats['chromosome'] = sumstats.chromosome.map(lambda x: variant_id_chr_map[x] if x in variant_id_chr_map else x)

    if set_variant_id:
        sumstats['rsid'] = sumstats.chromosome.map(lambda x: variant_id_chr_map[x] if x in variant_id_chr_map else x).str.cat(
            [sumstats.position.astype(str), sumstats.allele1, sumstats.allele2], sep=args.variant_id_sep)
        if grch38:
            sumstats['rsid'] = 'chr' + sumstats.rsid

    if maf_col in sumstats.columns:
        sumstats = sumstats.rename(index=str, columns={maf_col: 'maf'})
    elif freq_col in sumstats.columns:
        sumstats['maf'] = np.where(sumstats[freq_col] < 0.5, sumstats[freq_col], 1.0 - sumstats[freq_col])
    else:
        raise ValueError('sumstats should have maf/freq column.')

    # this happens when maf_col is misspecified (actually freq) or freq_col == 'maf'
    if np.any(sumstats.maf > 0.5):
        sumstats['maf'] = np.where(sumstats.maf < 0.5, sumstats.maf, 1.0 - sumstats.maf)
    if np.any(sumstats.maf <= 0):
        logger.warning('{} SNPs excluded due to MAF <= 0'.format(np.sum(sumstats.maf <= 0)))
        sumstats = sumstats.loc[sumstats.maf > 0, :]

    if flip_col in sumstats.columns:
        sumstats = sumstats.rename(index=str, columns={flip_col: 'flip'})
    else:
        sumstats['flip'] = 0
    if flip_beta:
        sumstats['beta'] = -sumstats.beta

    if p_col not in sumstats.columns:
        # assume standard linear regression
        sumstats['p'] = 2 * stats.norm.sf(np.abs(sumstats.beta / sumstats.se))
    elif p_col != 'p':
        sumstats = sumstats.rename(index=str, columns={p_col: 'p'})

    if scale_se_by_pval:
        se = np.abs(sumstats.beta / stats.norm.ppf(sumstats.p.astype(float) / 2))
        se[(sumstats.beta == 0) | ~np.isfinite(se)] = sumstats.se[(sumstats.beta == 0) | ~np.isfinite(se)]
        logger.info("{} SNPs are scaled (--scale-se-by-pval)".format(np.sum(~np.isclose(sumstats.se, se))))
        sumstats['se'] = se

    sumstats['chisq'] = (sumstats.beta / sumstats.se)**2
    sumstats = sumstats.dropna(subset=['beta', 'se', 'p'])
    return sumstats

def filter_sumstat(df, bed):
    for f in bed:
        df = df.loc[((df.chromosome == f.chrom) & (df.position >= f.start) & (df.position <= f.end)), :]
    return df


def generate_bed(sumstats,
                 p_threshold=5e-8,
                 maf_threshold=0,
                 window=0,
                 no_merge=False,
                 grch38=False,
                 exclude_MHC=False,
                 MHC_start=25e6,
                 MHC_end=34e6,
                 wdl=False,
                 min_p_threshold=None,
                 max_region_width=np.inf,
                 window_shrink_ratio=0.9):
    chisq_threshold = sp.stats.norm.ppf(p_threshold / 2) ** 2
    max_chisq_threshold = (sp.stats.norm.ppf(min_p_threshold / 2)**2) if min_p_threshold is not None else None

    df = pd.concat(map(lambda x: x[x.chisq > chisq_threshold], sumstats)).sort_values('chisq', ascending=False)
    bed_frames = []
    lead_snps = []

    build = 'hg19' if not grch38 else 'hg38'

    # exclude significant SNPs in the MHC region
    MHC_idx = (df.chromosome.map(CHROM_MAPPING_INT) == 6) & (df.position >= MHC_start) & (df.position <= MHC_end)
    if exclude_MHC and np.sum(MHC_idx) > 0:
        logger.warning('{} significant SNPs excluded due to --exclude-MHC'.format(np.sum(MHC_idx)))
        df = df.loc[~MHC_idx, :]
    # exclude by MAF
    if maf_threshold > 0:
        maf_idx = df.maf < maf_threshold
        logger.warning('{} significant SNPs excluded due to --maf-threshold'.format(np.sum(maf_idx)))
        df = df.loc[~maf_idx, :]
    if not wdl and len(df.index) == 0:
        raise RuntimeError('No signifcant SNPs found.')

    while len(df.index) > 0:
        lead_snp = df.iloc[0, :]
        chrom, start, end = lead_snp.chromosome, lead_snp.position - window, lead_snp.position + window
        chr_start, chr_end = pybedtools.chromsizes(build)['chr' + CHROM_MAPPING_STR[chrom]]
        start = max(start, chr_start)
        end = min(end, chr_end)
        df = df.loc[~((df.chromosome == chrom) & (df.position >= start) & (df.position <= end)), :]

        # Skip a lead SNP with p < min_p_threshold
        if max_chisq_threshold is not None and lead_snp.chisq > max_chisq_threshold:
            continue

        bed_frames.append(
            pd.DataFrame([[CHROM_MAPPING_INT[chrom], int(start), int(end)]], columns=['chrom', 'start', 'end']))
        lead_snps.append(lead_snp)

    if len(bed_frames) > 0:
        bed = BedTool.from_dataframe(pd.concat(bed_frames).sort_values(['chrom', 'start']))
        lead_snps = pd.concat(lead_snps, axis=1).T
        if not no_merge:
            print(len(sumstats[0].index))
            print("Window size {}".format(window))
            print("Merging regions. Before")
            print(bed)
            bed = bed.merge().saveas()
            bed_oversized = bed.filter(lambda x: len(x) > max_region_width).saveas()
            print("Oversized regions")
            print(bed_oversized)
            if (len(bed_oversized)) > 0:
                # keep in-size regions
                bed1 = bed.filter(lambda x: len(x) <= max_region_width).saveas()
                
                # filter to oversized regions
                sumstats = map(lambda x: filter_sumstat(x, bed_oversized), sumstats)

                # recursion with current window size * window_shrink_ratio
                bed2, lead_snps2 = generate_bed(
                    sumstats,
                    p_threshold=p_threshold,
                    maf_threshold=maf_threshold,
                    window=window * window_shrink_ratio,
                    no_merge=no_merge,
                    grch38=grch38,
                    exclude_MHC=exclude_MHC,
                    MHC_start=MHC_start,
                    MHC_end=MHC_end,
                    wdl=wdl,
                    min_p_threshold=min_p_threshold,
                    max_region_width=max_region_width,
                    window_shrink_ratio=window_shrink_ratio)

                bed = bed1.cat(bed2).saveas()
                lead_snps = pd.concat([lead_snps, lead_snps2], axis=0)
            print("After:")
            print(bed)
    else:
        bed = BedTool.from_dataframe(pd.DataFrame(columns=['chrom', 'start', 'end']))
        lead_snps = df

    return bed, lead_snps


def output_z(df, prefix, boundaries, grch38=False, no_output=True, extra_cols=None):
    i = int(df.name)
    outname = prefix + '.chr{}.{}-{}.z'.format(CHROM_MAPPING_STR[df.chromosome.iloc[0]], boundaries[i] % CHROM_CONSTANT,
                                               boundaries[i + 1] % CHROM_CONSTANT)
    if grch38:
        df['chromosome'] = 'chr' + df.chromosome

    output_cols = FINEMAP_COLUMNS
    if extra_cols is not None:
        output_cols = output_cols + extra_cols

    if not no_output:
        logger.info("Writing z file: " + outname)
        df[output_cols].to_csv(outname, sep=' ', float_format='%.6g', na_rep='NA', index=False)
    return outname


def output_tasks(args,
                 input_z,
                 region_size,
                 prefix,
                 gsdir,
                 localdir,
                 input_samples,
                 input_incl_samples,
                 n_samples,
                 var_y=None,
                 load_yty=False,
                 yty=None,
                 phi=None):

    if args.max_region_size < np.inf:
        print(region_size)
        idx = region_size < args.max_region_size
        input_z = input_z.loc[idx]
        logger.warning('{} regions excluded due to region_size <= {}'.format(np.sum(~idx), args.max_region_size))

    n_tasks = len(input_z)
    gsdir = pd.Series([gsdir] * n_tasks)
    localdir = pd.Series([localdir] * n_tasks)

    # ldstore
    input_base = input_z.str.slice(stop=-2)
    input_ld = input_base.str.cat(['.ld.bgz'] * n_tasks)
    out_bcor = input_base.str.cat(['.bcor'] * n_tasks)

    # finemap
    input_samples = [input_samples] * n_tasks
    input_incl_samples = [input_incl_samples] * n_tasks
    out_snp = input_base.str.cat(['.snp'] * n_tasks)
    out_config = input_base.str.cat(['.config'] * n_tasks)
    out_cred = input_base.str.cat(['.cred'] * n_tasks)
    out_log = input_base.str.cat(['.log_sss'] * n_tasks)
    n_samples = [n_samples] * n_tasks
    phi = [phi] * n_tasks

    # susie
    if load_yty:
        # gsdir should be added here to distinguish from *None*.
        input_yty = gsdir.str.cat(input_base.str.cat(['.yty'] * n_tasks).values)
    else:
        input_yty = [None] * n_tasks
    out_susie_snp = input_base.str.cat(['.susie.snp'] * n_tasks)
    out_susie_cred = input_base.str.cat(['.susie.cred'] * n_tasks)
    out_susie_log = input_base.str.cat(['.susie.log'] * n_tasks)
    out_susie_rds = input_base.str.cat(['.susie.rds'] * n_tasks)
    var_y = [var_y] * n_tasks
    yty = [yty] * n_tasks
    if args.dominant:
        dominant = [True] * n_tasks
    else:
        dominant = [None] * n_tasks

    ldstore_tasks = pd.DataFrame(
        OrderedDict((
            ('--input INPUT_Z', localdir.str.cat(input_z.values)),
            ('--input INPUT_SAMPLES', input_samples),
            ('--input INPUT_INCL_SAMPLES', input_incl_samples),
            ('--output OUT_BCOR', localdir.str.cat(out_bcor.values)),
            ('--output OUT_LD', localdir.str.cat(input_ld.values))
        )))
    finemap_tasks = pd.DataFrame(
        OrderedDict((
            ('--input INPUT_Z', gsdir.str.cat(input_z.values)),
            ('--input INPUT_LD', gsdir.str.cat(input_ld.values)),
            ('--output OUT_SNP', gsdir.str.cat(out_snp.values)),
            ('--output OUT_CONFIG', gsdir.str.cat(out_config.values)),
            ('--output OUT_CRED', gsdir.str.cat(out_cred.values)),
            ('--output OUT_LOG', gsdir.str.cat(out_log.values)),
            ('--env N_SAMPLES', n_samples),
            ('--env PHI', phi)
        )))
    susie_tasks = pd.DataFrame(
        OrderedDict((
            ('--input INPUT_Z', gsdir.str.cat(input_z.values)),
            ('--input INPUT_LD', gsdir.str.cat(input_ld.values)),
            ('--input INPUT_YTY', input_yty),
            ('--output OUT_SNP', gsdir.str.cat(out_susie_snp.values)),
            ('--output OUT_CRED', gsdir.str.cat(out_susie_cred.values)),
            ('--output OUT_LOG', gsdir.str.cat(out_susie_log.values)),
            ('--output OUT_RDS', gsdir.str.cat(out_susie_rds.values)),
            ('--env N_SAMPLES', n_samples),
            ('--env VAR_Y', var_y),
            ('--env YTY', yty),
            ('--env DOMINANT', dominant)
        )))

    if not args.no_ldstore:
        ldstore_tasks.to_csv(prefix + '.ldstore.tasks.txt', sep='\t', index=False)
    finemap_tasks.to_csv(prefix + '.finemap.tasks.txt', sep='\t', index=False)
    susie_tasks.to_csv(prefix + '.susie.tasks.txt', sep='\t', index=False)


def dsub(job_name,
         project,
         regions,
         machine_type,
         image,
         script,
         tasks,
         logging,
         preemptible=False,
         mount=None,
         env=None,
         after=None,
         dsub_shell_script=None,
         submit_jobs=False):
    cmd = [
        'dsub',
        '--provider', 'google-v2',
        '--project', project,
        '--regions', regions,
        '--machine-type', machine_type,
        '--image', image,
        '--name', job_name,
        '--script', script,
        '--tasks', tasks,
        '--logging', logging,
        '--disk-size', '100'
    ]

    if preemptible:
        cmd += ['--preemptible']
    if mount is not None:
        cmd += ['--mount', mount]
    if env is not None:
        if isinstance(env, str):
            env = [env]
        for e in env:
            cmd += ['--env', e]
    if after is not None:
        cmd += ['--after', after]

    logger.info("# " + job_name + " submission cmd:")
    formatted_cmd = ' '.join(cmd).replace(' --', ' \\\n--')
    logger.info(formatted_cmd)

    if dsub_shell_script is not None:
        with open(dsub_shell_script, 'w') as f:
            f.write('#!/bin/bash\n')
            f.write(formatted_cmd + '\n')
        make_executable(dsub_shell_script)

    if submit_jobs:
        if after is None:
            return run_command(cmd, stderr=None).strip()
        else:
            return run_command_bg(cmd)


def dsub_ldstore(args, prefix, gsdir, project, regions, bgen_bucket, bgen_dirname, bgen_fname_format,
                 submit_jobs=False):
    return dsub(
        'ldstore-' + prefix,
        project=project,
        regions=regions,
        machine_type=args.ldstore_machine_type,
        image=args.ldstore_image,
        script=args.ldstore_script,
        tasks=prefix + '.ldstore.tasks.txt',
        dsub_shell_script=prefix + '.ldstore.dsub.sh',
        logging=gsdir + 'logging',
        preemptible=args.preemptible,
        mount=bgen_bucket,
        env=['BGEN_DIRNAME='+bgen_dirname, 'BGEN_FNAME_FORMAT='+bgen_fname_format],
        submit_jobs=args.submit_jobs)


def dsub_finemap(args, prefix, gsdir, project, regions, submit_jobs=False, after=None):
    return dsub(
        'finemap-' + prefix,
        project=project,
        regions=regions,
        machine_type=args.finemap_machine_type,
        image=args.finemap_image,
        script=args.finemap_script,
        tasks=prefix + '.finemap.tasks.txt',
        dsub_shell_script=prefix + '.finemap.dsub.sh',
        logging=gsdir + 'logging',
        preemptible=args.preemptible,
        after=after,
        submit_jobs=args.submit_jobs)


def dsub_susie(args, prefix, gsdir, project, regions, submit_jobs=False, after=None):
    return dsub(
        'susie-' + prefix,
        project=project,
        regions=regions,
        machine_type=args.susie_machine_type,
        image=args.susie_image,
        script=args.susie_script,
        tasks=prefix + '.susie.tasks.txt',
        dsub_shell_script=prefix + '.susie.dsub.sh',
        logging=gsdir + 'logging',
        preemptible=args.preemptible,
        after=after,
        submit_jobs=args.submit_jobs)


def main(args):
    # remove x chromosome from the key unless --x-chromosome is specified
    if not args.x_chromosome:
        max_chrom_int = 22
        for k in ['X', '23']:
            del CHROM_MAPPING_INT[k]
            del CHROM_MAPPING_STR[k]
    else:
        max_chrom_int = 23

    # read sumstats
    var_id_chr_map = {}

    if args.set_variant_id_map_chr:
        var_id_chr_map = { m[0].strip():m[1].strip() for m in map(lambda x: x.split("="), args.set_variant_id_map_chr.split(",")) }

    sumstats = map(
        lambda (i, x): read_sumstats(
            x,
            rsid_col=args.rsid_col[i],
            chromosome_col=args.chromosome_col[i],
            position_col=args.position_col[i],
            allele1_col=args.allele1_col[i],
            allele2_col=args.allele2_col[i],
            maf_col=args.maf_col[i],
            freq_col=args.freq_col[i],
            beta_col=args.beta_col[i],
            se_col=args.se_col[i],
            flip_col=args.flip_col[i],
            p_col=args.p_col[i],
            delimiter=args.delimiter[i],
            recontig=args.recontig[i],
            set_variant_id=args.set_variant_id[i],
            variant_id_chr_map=var_id_chr_map,
            flip_beta=args.flip_beta[i],
            grch38=args.grch38,
            scale_se_by_pval=args.scale_se_by_pval[i],
            extra_cols=args.extra_cols
        ), enumerate(args.sumstats))


    if args.bed is None:
        logger.info('Generating bed')
        merged_bed, lead_snps = generate_bed(sumstats, p_threshold=args.p_threshold, maf_threshold=args.maf_threshold,
                                            window=args.window, no_merge=args.no_merge,
                                            grch38=args.grch38, exclude_MHC=args.exclude_MHC,
                                            MHC_start=args.MHC_start, MHC_end=args.MHC_end, wdl=args.wdl,
                                            min_p_threshold=args.min_p_threshold,
                                            max_region_width=args.max_region_width,
                                            window_shrink_ratio=args.window_shrink_ratio)
        lead_snps.to_csv(args.out + '.lead_snps.txt', sep='\t', index=False)
    else:
        i = 0
        logger.info('Loading user-supplied bed: ' + args.bed[i])
        bed = pd.read_csv(args.bed[i],
                          delim_whitespace=True,
                          header=None,
                          names=['chromosome', 'start', 'end'],
                          dtype={'chromosome': str})
        bed['chromosome'] = bed.chromosome.map(CHROM_MAPPING_INT)
        merged_bed = BedTool.from_dataframe(bed)
        if not args.no_merge:
            merged_bed = merged_bed.merge()
            merged_bed = BedTool.from_dataframe(bed)

    logger.info(merged_bed)

    # write had results indicator file for WDL purposes
    if args.wdl:
        with open(args.out + "_had_results", 'w') as o:
            if (merged_bed.count() == 0):
                logger.info("No significant regions identified.")
                o.write("False\n")
                return
            else:
                o.write("True\n")

    build = 'hg19' if not args.grch38 else 'hg38'
    chromsizes = pd.DataFrame.from_dict(
        pybedtools.chromsizes(build), orient='index',
        columns=['start', 'end']).loc[['chr' + CHROM_MAPPING_STR[str(i)] for i in range(1, max_chrom_int + 1)], :]
    chromsizes['chromosome'] = range(1, max_chrom_int+1)
    chromsizes = chromsizes[['chromosome', 'start', 'end']]

    unique_chroms = merged_bed.to_dataframe().iloc[0, :].unique()
    all_bed = BedTool.from_dataframe(chromsizes[chromsizes.chromosome.isin(unique_chroms)])
    all_bed = all_bed.subtract(merged_bed).cat(
        merged_bed, postmerge=False).to_dataframe().sort_values(['chrom', 'start'])
    logger.info(all_bed)
    merged_bed = merged_bed.to_dataframe()

    boundaries = pd.concat([convert_chrpos(all_bed.chrom, all_bed.start),
                            convert_chrpos(all_bed.chrom, all_bed.end)]).sort_values().unique()
    logger.debug(boundaries)

    if not args.no_merge:
        # assign regions
        sumstats = map(
            lambda x: x.assign(region=pd.cut(x.chrpos, boundaries, right=False, labels=False, include_lowest=True)),
            sumstats)
        if not args.null_region:
            sig_regions = pd.cut(convert_chrpos(merged_bed.chrom, merged_bed.start),
                                 boundaries,
                                 right=False,
                                 labels=False,
                                 include_lowest=True).unique()
            sumstats = map(lambda x: x[x.region.isin(sig_regions)], sumstats)
        else:
            merged_bed = all_bed

    if args.bed is None:
        merged_bed.to_csv(args.out + '.bed', sep='\t', index=False, header=False)

    for i, x in enumerate(sumstats):
        if args.extract is not None:
            if 'v' not in x.columns:
                x['v'] = x.chromosome.map(CHROM_MAPPING_STR) + ':' + x.position.astype(
                    str) + ':' + x.allele1 + ':' + x.allele2
            snplist = pd.read_csv(args.extract,
                                  delim_whitespace=True,
                                  header=None,
                                  compression='gzip' if args.extract.endswith('gz') else 'infer').iloc[:, 0]
            logger.warning('{} SNPs excluded due to --extract'.format(np.sum(~x.v.isin(snplist))))
            x = x.loc[x.v.isin(snplist), :]

        if not args.no_merge:
            input_z = x.groupby('region').apply(output_z,
                                                prefix=args.prefix[i],
                                                boundaries=boundaries,
                                                grch38=args.grch38,
                                                no_output=args.no_output,
                                                extra_cols=args.extra_cols)

            region_size = x.groupby('region').count().rsid
        else:
            # boundaries are ordered wrongly if there are overlapping regions
            boundaries = pd.concat(
                [convert_chrpos(merged_bed.chrom, merged_bed.start),
                 convert_chrpos(merged_bed.chrom, merged_bed.end)],
                axis=1).values.flatten()
            logger.debug(boundaries)

            # dirty hack for now -- pd.cut doesn't accept variants assigned to multiple regions
            input_z = []
            region_size = []
            for j, row in merged_bed.iterrows():
                start = convert_chrpos(row.chrom, row.start)
                end = convert_chrpos(row.chrom, row.end)
                idx = (start <= x.chrpos) & (x.chrpos <= end)
                # region values should be consistent with pd.cut outputs
                xx = x.loc[idx, :]
                xx.name = 2 * j
                input_z.append(
                    output_z(xx,
                             prefix=args.prefix[i],
                             boundaries=boundaries,
                             grch38=args.grch38,
                             no_output=args.no_output,
                             extra_cols=args.extra_cols))
                region_size.append(np.sum(idx))

            input_z = pd.Series(input_z)
            region_size = pd.Series(region_size)

        if not args.no_upload:
            logger.info("Uploading z files")
            if args.bed is None:
                bed_and_leadsnps = [args.out + '.bed', args.out + '.lead_snps.txt']
            else:
                bed_and_leadsnps = [args.bed[0]]
            run_command(['gsutil', '-m', 'cp'] + bed_and_leadsnps + input_z.values.tolist() + [args.gsdir[i]])

        if args.wdl:
            logger.info("--wdl is specified. No task/dsub files were output.")
            return

        logger.info("Writing task files")
        output_tasks(args, input_z, region_size, args.prefix[i], args.gsdir[i], args.localdir[i], args.input_samples[i],
                     args.input_incl_samples[i], args.n_samples[i], args.var_y[i], args.load_yty[i], args.yty[i],
                     args.phi[i])

        if args.localdir[i].startswith('gs://'):
            if not args.no_ldstore:
                job_id = dsub_ldstore(args, args.prefix[i], args.gsdir[i], args.project[i], args.regions[i],
                                      args.bgen_bucket[i], args.bgen_dirname[i], args.bgen_fname_format[i],
                                      args.submit_jobs)
            else:
                job_id = None
            dsub_finemap(
                args, args.prefix[i], args.gsdir[i], args.project[i], args.regions[i], args.submit_jobs, after=job_id)
            dsub_susie(
                args, args.prefix[i], args.gsdir[i], args.project[i], args.regions[i], args.submit_jobs, after=job_id)
        else:
            dsub_finemap(args, args.prefix[i], args.gsdir[i], args.project[i], args.regions[i], submit_jobs=False)
            dsub_susie(args, args.prefix[i], args.gsdir[i], args.project[i], args.regions[i], submit_jobs=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=str)
    parser.add_argument('--p-threshold', type=float, default=5e-8)
    parser.add_argument(
        '--min-p-threshold',
        type=float,
        help=(
            "Minimum lead variant p-value for inclusion in finemapping. "
            "Useful for adding less significant regions after genome-wide significant finemapping has already been done"
        ))
    parser.add_argument('--maf-threshold', type=float, default=0, help='MAF threshold for lead variants')
    parser.add_argument('--max-region-size',
                        type=int,
                        default=np.inf,
                        help='Maximum number of variants limit for a region to finemap (only applies to task files)')
    parser.add_argument('--window', type=int, default=1.5e6)
    parser.add_argument('--max-region-width',
                        type=int,
                        default=np.inf,
                        help='Maximum width of finemap regions after possible merge')
    parser.add_argument('--window-shrink-ratio',
                        type=float,
                        default=0.9,
                        help='Ratio to recursively shrink a flanking window size')
    parser.add_argument('--no-merge', action='store_true', help='Do not merge overlapped regions')
    parser.add_argument('--null-region', action='store_true')
    parser.add_argument('--no-upload', action='store_true')
    parser.add_argument('--no-output', action='store_true')
    parser.add_argument('--no-ldstore', action='store_true')
    parser.add_argument('--submit-jobs', action='store_true')
    parser.add_argument('--exclude-MHC', action='store_true')
    parser.add_argument('--MHC-start', type=int, default=25e6)
    parser.add_argument('--MHC-end', type=int, default=34e6)
    parser.add_argument('--extract', type=str, help='Extract specific set of variants to fine-map')
    parser.add_argument('--dominant', action='store_true')
    parser.add_argument('--x-chromosome', action='store_true')
    parser.add_argument('--variant-id-sep',
                        type=str,
                        default='_',
                        help='Separator for variant ID when --set-variant-id is specified.')
    parser.add_argument('--wdl', action='store_true', help='Output extra files for wdl')

    # sumstats settings
    parser.add_argument('--json', type=str, nargs='+')
    parser.add_argument('--sumstats', type=str, nargs='+')
    parser.add_argument('--prefix', type=str, nargs='+')
    parser.add_argument('--bed', type=str, nargs='+')
    parser.add_argument('--rsid-col', type=str, default='rsid')
    parser.add_argument('--chromosome-col', type=str, default='chromosome')
    parser.add_argument('--position-col', type=str, default='position')
    parser.add_argument('--allele1-col', type=str, default='allele1')
    parser.add_argument('--allele2-col', type=str, default='allele2')
    parser.add_argument('--maf-col', type=str, default='maf')
    parser.add_argument('--freq-col', type=str, default='freq')
    parser.add_argument('--beta-col', type=str, default='beta')
    parser.add_argument('--se-col', type=str, default='se')
    parser.add_argument('--flip-col', type=str, default='flip')
    parser.add_argument('--p-col', '-p', type=str, default='p')
    parser.add_argument('--delimiter', type=str, default='\s+', help='Delimiter of sumstats')
    parser.add_argument('--extra-cols', type=str, nargs='+', help='Extra columns to output in .z files.')
    parser.add_argument('--recontig', action='store_true', default=False)
    parser.add_argument('--set-variant-id', action='store_true', default=False)
    parser.add_argument('--set-variant-id-map-chr', type=str,
                        help="Comma separated list of chromosome id remapping to be done for the generated variant. "
                        "E.g. to map 23 to X and 24 to MT specify 23=X,24=MT."
                        "This is useful if there is a mismatch in variant id in bgen and summary stats."
                        )
    parser.add_argument('--flip-beta', action='store_true', default=False)
    parser.add_argument('--grch38', action='store_true', default=False)
    parser.add_argument(
        '--scale-se-by-pval',
        action='store_true',
        default=False,
        help='Get Z score from p-value instead of using beta and se. '
        'Should be used when SPA approximation has been used to generate p-value (e.g. SAIGE)')

    JSON_PARAMS = [
        'rsid_col', 'chromosome_col', 'position_col', 'allele1_col', 'allele2_col', 'maf_col', 'freq_col', 'beta_col',
        'se_col', 'flip_col', 'p_col', 'delimiter', 'recontig', 'set_variant_id', 'flip_beta', 'scale_se_by_pval', 'project',
        'regions', 'bgen_bucket', 'bgen_dirname', 'bgen_fname_format', 'var_y', 'load_yty', 'yty', 'phi'
    ]

    # task parameters (usually set by json)
    parser.add_argument('--gsdir', type=str, nargs='+')
    parser.add_argument('--localdir', type=str, nargs='+')
    parser.add_argument('--input-samples', type=str, nargs='+')
    parser.add_argument('--input-incl-samples', type=str, nargs='+')
    parser.add_argument('--n-samples', '-n', type=int, nargs='+')
    parser.add_argument('--var-y', type=float, nargs='+', default=None)
    parser.add_argument('--load-yty', action='store_true', default=False)
    parser.add_argument('--yty', type=float, nargs='+', default=None)
    parser.add_argument('--phi', type=float, nargs='+', default=None)

    # dsub settings
    parser.add_argument('--project', type=str, default='encode-uk-biobank-restrict', nargs='+')
    parser.add_argument('--regions', type=str, default='us-central1', nargs='+')
    parser.add_argument('--preemptible', action='store_true')
    parser.add_argument('--bgen-bucket',
                        type=str,
                        default='BGEN_BUCKET=gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183',
                        nargs='+')
    parser.add_argument('--bgen-dirname', type=str, default='imputed', nargs='+')
    parser.add_argument('--bgen-fname-format', type=str, default='ukb_imp_chr{}_v3.bgen', nargs='+')
    parser.add_argument('--ldstore-machine-type', type=str, default='n1-standard-16')
    parser.add_argument('--ldstore-image', type=str, default='gcr.io/encode-uk-biobank-restrict/finemap-suite:0.6')
    parser.add_argument(
        '--ldstore-script',
        type=str,
        default='/humgen/atgu1/fs03/mkanai/workspace/201901_XFinemap/xfinemap/docker/ldstore/dsub_ldstore.py')
    parser.add_argument('--finemap-machine-type', type=str, default='n1-highmem-4')
    parser.add_argument('--finemap-image', type=str, default='gcr.io/encode-uk-biobank-restrict/finemap:1.3.1')
    parser.add_argument(
        '--finemap-script',
        type=str,
        default='/humgen/atgu1/fs03/mkanai/workspace/201901_XFinemap/xfinemap/docker/finemap/dsub_finemap.py')
    parser.add_argument('--susie-machine-type', type=str, default='n1-highmem-16')
    parser.add_argument('--susie-image',
                        type=str,
                        default='gcr.io/encode-uk-biobank-restrict/susie:0.8.1.0521.save-rds.yty.dominant')
    parser.add_argument(
        '--susie-script',
        type=str,
        default='/humgen/atgu1/fs03/mkanai/workspace/201901_XFinemap/xfinemap/docker/susie/dsub_susie.py')

    args = parser.parse_args()

    n_sumstats = len(args.sumstats)

    args_dict = vars(args)
    defaults = vars(parser.parse_args(''))

    if args.json is not None:
        if n_sumstats != len(args.json):
            raise ValueError("--sumstats and --json should have the same length.")

        json_dict = defaultdict(list)

        def update_json_dict(file):
            d = json.load(open(file))
            if 'localdir' not in d and 'gsdir' in d:
                d.update({'localdir': d['gsdir']})
            for k, v in d.items():
                json_dict[k].append(v)
            for k in np.setdiff1d(JSON_PARAMS, d.keys()):
                json_dict[k].append(defaults[k])

        map(lambda x: update_json_dict(x), args.json)
        args_dict.update(json_dict)

    if args.prefix is None:
        args.prefix = map(
            lambda x: os.path.splitext(os.path.basename(x if not x.endswith('gz') else os.path.splitext(x)[0]))[0],
            args.sumstats)
    if args.gsdir is not None:
        args.gsdir = map(lambda x: x + '/' if not x.endswith('/') else x, args.gsdir)
    if args.localdir is not None:
        args.localdir = map(lambda x: x + '/' if not x.endswith('/') else x, args.localdir)

    if args.bed is not None and not isinstance(args.bed, list):
        args.bed = [args.bed] * n_sumstats

    args_dict.update({
        k: [v] * n_sumstats if not (isinstance(v, list)) else v
        for k, v in args_dict.items()
        if k in JSON_PARAMS
    })

    if args.delimiter is not None:
        def update_delimiter(delimiter):
            STANDARD_DELIMITERS = ['\s', '\s+', '\t', ' ']
            MAGIC_WORDS = ['SINGLE_WHITESPACE', 'WHITESPACE', 'TAB', 'SPACE']
            if delimiter == 'SINGLE_WHITESPACE':
                return '\s'
            elif delimiter == 'WHITESPACE':
                return '\s+'
            elif delimiter == 'TAB':
                return '\t'
            elif delimiter == 'SPACE':
                return ' '
            elif delimiter not in STANDARD_DELIMITERS:
                logger.warning(
                    '--delimiter %s does not seem a standard delimiter nor match any of magic words (%s).'.format(
                        args.delimiter, ','.join(MAGIC_WORDS)))
                return delimiter
        args.delimiter = map(lambda x: update_delimiter(x), args.delimiter)
    else:
        raise ValueError('--delimiter should be specified.')

    if not args.no_upload:
        len_check_params = [args.gsdir]
        if not args.wdl:
            len_check_params += [args.localdir, args.input_samples, args.input_incl_samples, args.n_samples, args.var_y]
        if args.gsdir is None:
            raise ValueError("--gsdir should be specified.")
        if not np.all(n_sumstats == np.array(map(len, len_check_params))):
            raise ValueError("Different length.")

    if args.null_region and args.no_merge:
        raise ValueError("--null-region and --no-merge cannot be specified at the same time.")

    if args.out is None:
        if n_sumstats == 1:
            args.out = args.prefix[0]
        else:
            raise ValueError("--out should be specified when # sumstats > 1.")

    if args.submit_jobs and 'GOOGLE_APPLICATION_CREDENTIALS' not in os.environ:
        raise RuntimeError("GOOGLE_APPLICATION_CREDENTIALS should be set in the environment variables.")

    # logging
    fhandler = FileHandler(args.out + '.log', 'w+')
    fhandler.setFormatter(log_format)
    logger.addHandler(fhandler)

    # https://github.com/bulik/ldsc/blob/master/ldsc.py#L589
    defaults = vars(parser.parse_args(''))
    non_defaults = [x for x in args_dict.keys() if args_dict[x] != defaults[x]]
    # non_defaults = [x for x in args_dict.keys()]
    call = "Call: \n"
    call += './{}.py \\\n'.format(os.path.basename(__file__))
    options = ['--' + x.replace('_', '-') + ' ' + str(args_dict[x]) + ' \\' for x in non_defaults]
    call += '\n'.join(options) #.replace('True', '').replace('False', '')
    call = call[0:-1] + '\n'
    logger.info(call)

    main(args)
