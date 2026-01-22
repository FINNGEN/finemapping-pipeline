# coding: utf-8
import argparse
import numpy as np
import hail as hl
import re
from hail.linalg import BlockMatrix
from hail.utils import new_temp_file
# from ukbb_pan_ancestry.resources.ld import get_ld_matrix_path, get_ld_variant_index_path

bucket = 'gs://pan-ukbb-refinery-ld'


def get_ld_matrix_path(pop: str, extension: str = 'bm'):
    if pop == 'FG':
        return f'gs://r5_data/ld/FG.r5.ldadj.{extension}'
    return f'{bucket}/ld_release/UKBB.{pop}.ldadj.{extension}'


def get_ld_variant_index_path(pop: str, extension: str = 'ht'):
    if pop == 'FG':
        return f'gs://r5_data/ld/FG.r5.ldadj.variant.{extension}'
    return f'{bucket}/ld_release/UKBB.{pop}.ldadj.variant.liftover.{extension}'


def get_diag_mat(diag_vec: BlockMatrix):
    x = diag_vec.T.to_numpy()
    diag_mat = np.identity(len(x)) * np.outer(np.ones(len(x)), x)
    return BlockMatrix.from_numpy(diag_mat)


def densify_mat(bm: BlockMatrix, sparsity):
    if sparsity == 'triangular':
        bm = bm.sparsify_triangle()
    elif sparsity == 'full':
        # re-densify
        bm = bm + bm.T - get_diag_mat(bm.diagonal())

    return bm


def export_ld(bm: BlockMatrix, out, delimiter=' '):
    bm_tmp = new_temp_file()
    bm = bm.checkpoint(bm_tmp, force_row_major=True)
    bm.export(
        bm_tmp,
        out,
        delimiter=delimiter,
    )

    return bm


def align_indices(value, idx, n_variants):
    mat = np.zeros((n_variants, n_variants))
    mat[np.ix_(idx, idx)] = value
    bm = BlockMatrix.from_numpy(mat)
    bm = bm.checkpoint(new_temp_file(), force_row_major=True)

    return bm


def main(args):

    hl.init(tmp_dir=args.tmpdir, local_tmpdir=args.tmpdir, default_reference='GRCh38')

    ht_z = hl.import_table(args.zfile_meta, delimiter=' ', impute=True)
    ht_z = ht_z.key_by(**hl.parse_variant(ht_z.rsid))

    hts = []
    for pop in args.pops:
        print(pop)
        # import zfile and join with index
        ht = hl.read_table(get_ld_variant_index_path(pop))
        ht = ht.annotate_globals(pop=pop)
        ht = ht_z.join(ht, 'inner')
        ht = ht.checkpoint(new_temp_file())
        hts.append(ht.select('idx').select_globals('pop', 'n_samples'))

        print(f'{ht.count()} variants found.')
        if ht.aggregate(hl.agg.any(hl.is_missing(ht.idx))):
            ht = ht.filter(hl.is_missing(ht.idx))
            ht.show()
            raise RuntimeError(f'{ht.count()} Indices are missing')

    # generate meta index
    ht = hl.Table.multi_way_zip_join(hts, 'row_field_name', 'global_field_name')
    ht = ht.transmute(idx=hl.map(lambda x: x.idx, ht.row_field_name))
    ht = ht.transmute_globals(pop_idx=hl.dict(
        hl.enumerate(hl.map(lambda x: x.pop, ht.global_field_name), index_first=False)),
                              n_samples=hl.map(lambda x: x.n_samples, ht.global_field_name))
    ht = ht.add_index(name='idx_meta')
    ht = ht.checkpoint(new_temp_file())

    if args.n_samples is None:
        pop_indices = [ht.pop_idx[pop] for pop in args.pops]
        args.n_samples = ht.n_samples[pop_indices].collect()

    n_variants = ht.count()
    numer = None
    denom = 0

    for pop, n_samples in zip(args.pops, args.n_samples):
        print(f'{n_samples} samples found.')
        bm = BlockMatrix.read(get_ld_matrix_path(pop))

        pop_idx = ht.pop_idx[pop].collect()[0]
        # keep defined idx for a pop
        ht_pop = ht.filter(hl.is_defined(ht.idx[pop_idx]))
        idx = ht_pop.idx[pop_idx].collect()
        idx2 = sorted(idx)
        print(np.all(np.diff(idx) > 0))
        print(np.all(np.diff(idx) >= 0))

        bm = bm.filter(idx2, idx2)
        if not np.all(np.diff(idx) > 0):
            order = np.argsort(idx)
            rank = np.empty_like(order)
            rank[order] = np.arange(len(idx))
            mat = bm.to_numpy()[np.ix_(rank, rank)]
            bm = BlockMatrix.from_numpy(mat)
        bm = densify_mat(bm, sparsity=args.sparsity)
        # bm = export_ld(bm, args.out_pattern.format(POP=pop), delimiter=args.delimiter)

        if args.meta:
            idx_meta = ht_pop.idx_meta.collect()
            bm = n_samples * align_indices(bm.to_numpy(), idx_meta, n_variants)

            numer = numer + bm if numer is not None else bm
            denom += n_samples

    if args.meta:
        # add a min float to avoid zero devision
        bm_meta = numer / denom

        # merge with meta zfile
        print("Meta")
        print(f'{ht_z.count()} variants found.')
        ht_z = ht_z.join(ht.select('idx_meta'), 'inner')
        ht_z.key_by().drop('locus', 'alleles', 'idx_meta').export(f'{args.out}.z', delimiter=' ')

        print(f'{ht_z.count()} variants remained.')
        idx_meta = ht_z.idx_meta.collect()
        bm_meta = bm_meta.filter(idx_meta, idx_meta)

        if args.sparsity == 'triangular':
            bm_meta = densify_mat(bm_meta, sparsity=args.sparsity)

        export_ld(bm_meta, f'{args.out}.ld.bgz', delimiter=args.delimiter)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pops', type=str, nargs='+', required=True, help='Population to export a LD matrix')
    parser.add_argument('--n-samples', type=int, nargs='+', help='Number of samples for each population')
    parser.add_argument('--zfile-meta',
                        type=str,
                        required=True,
                        help='Input z file from meta-analysis')
    parser.add_argument('--out',
                        type=str,
                        help='Output prefix')
    parser.add_argument('--delimiter', type=str, default=' ', help='Delimiter for output ld matrix')
    parser.add_argument('--sparsity',
                        type=str,
                        default='keep',
                        const='keep',
                        nargs='?',
                        choices=['keep', 'triangular', 'full'])
    parser.add_argument('--meta', action='store_true', help='Write weighted average of LD matrices for meta')
    parser.add_argument('--tmpdir', type=str, help='Tmporary directory for hail')

    args = parser.parse_args()

    if (args.n_samples is not None) and (len(args.pops) != len(args.n_samples)):
        raise ValueError('--n-samples should have the same length as --pops.')

    main(args)
