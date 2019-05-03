import atexit
import os
import os.path
import shutil
import stat
import sys
import tempfile
import numpy as np
import pandas as pd
import subprocess as sp


# https://github.com/hail-is/hail/blob/master/hail/python/hail/utils/misc.py
def new_local_temp_dir(suffix='', prefix='tmp', dir=None):
    local_temp_dir = tempfile.mkdtemp(suffix, prefix, dir)
    atexit.register(shutil.rmtree, local_temp_dir)
    return local_temp_dir


def run_command(args, stderr=sp.STDOUT):
    try:
        if sys.version_info < (3,):
            output = sp.check_output(args, stderr=stderr)
        else:
            output = sp.check_output(args, encoding='utf-8', stderr=stderr)
        print(output)
        return(output)
    except sp.CalledProcessError as e:
        print(e.output)
        raise e


def run_command_bg(args):
    try:
        return sp.Popen(args)
    except sp.CalledProcessError as e:
        print(e.output)
        raise e


def gzip(file, decompress=False, keep=False, bgzip=False, n_threads=None):
    if bgzip:
        gzip, ext = 'bgzip', '.bgz'
    else:
        gzip, ext = 'gzip', '.gz'

    cmd = [gzip, file]
    if decompress:
        outname, suffix = os.path.splitext(file)
        cmd += ['-d', '-S', suffix]
        mode = 'w'
    else:
        outname = file + ext
        mode = 'wb'

    if keep:
        cmd += ['-c']
        f = open(outname, mode)
        stdout = f
    else:
        stdout = None
    if bgzip and n_threads is not None:
        cmd += ['-@', n_threads]

    try:
        sp.check_call(cmd, stdout=stdout)
    except sp.CalledProcessError as e:
        print(e.output)
        raise e

    if keep:
        f.close()
    if bgzip and not keep and not decompress:
        os.rename(file + '.gz', outname)

    return outname


def gunzip(file, keep=False, bgzip=False, n_threads=None):
    return gzip(file, decompress=True, keep=keep, bgzip=bgzip, n_threads=n_threads)


# https://stackoverflow.com/questions/12791997/how-do-you-do-a-simple-chmod-x-from-within-python
def make_executable(path):
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2  # copy R bits to X
    os.chmod(path, mode)


def pad_arrays(arr1, arr2, pad_width1, pad_width2, constant_values):
    arr1 = np.pad(arr1, (0, pad_width1), 'constant', constant_values=constant_values)
    arr2 = np.pad(arr2, (0, pad_width2), 'constant', constant_values=constant_values)
    return (arr1, arr2)


def pad_2d_array(arr, pad_width1, pad_width2, constant_values):
    arr = np.pad(arr, ((0, pad_width1), (0, pad_width2)), 'constant', constant_values=constant_values)
    return (arr)


def pad_str_arrays(arr1, arr2, pad_width1, pad_width2):
    size1, size2 = len(arr1) + pad_width1, len(arr2) + pad_width2
    new_arr1 = [''] * size1
    new_arr1[:len(arr1)] = arr1.tolist()
    new_arr2 = [''] * size2
    new_arr2[:len(arr2)] = arr2.tolist()
    return (new_arr1, new_arr2)


def rearrange_array(arr, new_size, ix, default_value=np.nan):
    new_arr = np.tile(default_value, new_size)
    new_arr[ix] = arr
    return (new_arr)


def _varid_to_rsid(self, result1, result2):
    ix1 = pd.Series(result1.snp).str.startswith('rs').values
    ix2 = pd.Series(result2.snp).str.startswith('rs').values
    d1 = pd.Series(result1.snp[ix1], index=result1.varid[ix1]).to_dict()
    d2 = pd.Series(result2.snp[ix2], index=result2.varid[ix2]).to_dict()

    # non rsids mapped to varid
    ix1, ix2 = np.logical_not(ix1), np.logical_not(ix2)
    nd1 = pd.Series(result1.varid[ix1], index=result1.snp[ix1]).to_dict()
    nd2 = pd.Series(result2.varid[ix2], index=result2.snp[ix2]).to_dict()
    self._rsid_dict = dict(nd1.items() + nd2.items() + d1.items() + d2.items())
    snp = pd.Series(self.varid)
    snp = snp.map(self._rsid_dict).fillna(snp).values
    return (snp)