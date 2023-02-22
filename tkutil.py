import scipy.sparse
import pickle
import gzip
import pandas as pd
import numpy as np
import scipy.io
import os, sys, re
import logging
import argparse
import io
import numba
# __loggers = {}

@numba.njit('i4(f4[:],i4[:,:])')
def _search_sections_notnan(values, sections):
    index = 0
    start = end = -1
    for i, v in enumerate(values):
        if np.isnan(v) is False:
            if start < 0:
                start = i
        elif start >= 0:
            end = i
            sections[:,index] = [start, end]
            start = -1
            index += 1
    if start >= 0:
        sections[:,index] = [start, len(values)]
        index += 1
    return index

@numba.njit('i4(f4[:],i4[:,:])')
def _search_sections_notzero(values, sections):
    index = 0
    start = end = -1
    for i, v in enumerate(values):
        if np.isnan(v) is False and v > 0:
            if start < 0:
                start = i
        elif start >= 0:
            end = i
            sections[:,index] = [start, end]
            start = -1
            index += 1
    if start >= 0:
        sections[:,index] = [start, len(values)]
        index += 1
    return index                

def smoothen_data_with_null(values, method='savgol', abovezero=False, **kwargs)->np.array:
    """Smoothen data with NaN or zero values.
    Args:
        values (1d np.array): Data undergone smoothing.
        method (str, optional): Smoothing algorithm. Defaults to 'savgol'.
        window : smoothing width length
        
        savgol_filter
        mode : 'mirror', 'nearest', 'constant', 'wrap'
        window_length : = window
        polyorder : less than window_length
    """
    import numba
    import scipy.signal
    # import scipy.interpolate
    # import statsmodels.nonparametric.smoothers_lowess

    window = kwargs.get('window', 11)
    polyorder = kwargs.get('polyorder', 3)
    remove_zero = abovezero
    values = np.array(values, dtype=np.float32)
    
    if method == 'savgol':
        if window % 2 == 0:
            window += 1
        if polyorder < 3: polyorder = 3
    else:
        raise Exception(f'not supported smoothing method {method}')

    sections = np.full((2, values.size // 2 + 1), -1, dtype=np.int32)
    if remove_zero:    
        n_sections = _search_sections_notzero(values, sections)
    else:
        n_sections = _search_sections_notnan(values, sections)
        
    smoothened = np.full(values.size, np.nan, dtype=np.float32)        
    for i in range(n_sections):
        sec = sections[:,i].reshape(-1)#for i, sec in enumerate(sections):
        if sec[0] < 0:
            break
        # print(i, sec)
        subval = values[sec[0]:sec[1]]
        if window > 0 and subval.size >= window:
            if method == 'savgol':
                subval = scipy.signal.savgol_filter(subval, window, polyorder)
        smoothened[sec[0]:sec[1]] = subval
        pass
    return smoothened

def check_file(filename, date_stat=None):
    """Examine processed file. If the file exists newer than other file, return True
"""
    flag_exist =  os.path.isfile(filename) and os.path.getsize(filename) > 1024
    if not flag_exist:
        return False
    if date_stat is None:
        return True
    if isinstance(date_stat, str) and os.path.isfile(date_stat):
        date_stat = int(os.stat(date_stat).st_mtime)
    mtime = int(os.stat(filename).st_mtime)
    if isinstance(date_stat, int):
        return mtime >= date_stat
    if hasattr(date_stat, 'st_mtime'):
        return mtime >= int(date_stat.st_mtime)
    return True

def get_logger(name=None, stdout=True, logfile=None):
    if name is None:
        name = sys._getframe().f_code.co_name
        pass
    
    logger = logging.getLogger(name)
    # set logging
    for h in logger.handlers:
        h.removeHandler(h)
    def _set_log_handler(logger, handler):#, verbose):
        handler.setFormatter(logging.Formatter('%(asctime)s %(name)s:%(lineno)s %(funcName)s [%(levelname)s]: %(message)s'))
        logger.addHandler(handler)
        return logger
    if logfile is not None:
        _set_log_handler(logger, logging.FileHandler(logfile))
    else:
        stdout = True
    if stdout:
        _set_log_handler(logger, logging.StreamHandler())
    # logger.setLevel(logging.ERROR)
    logger.propagate = False
    return logger

def instanciate_standard_logger(name=None):
    if name is None: name = __name__
    return get_logger(name, stdout=True, logfile='.run.log')

def read_chromosomes(sambamfile, logger=None):
    import collections
    cmd = 'samtools', 'view', '-H', sambamfile
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    chromosomes = collections.OrderedDict()
    with io.TextIOWrapper(proc.stdout, encoding='ascii') as fi:
        for line in fi:
            m = re.match('@SQ\\s+SN:(\\S+)\\s+LN:(\\d+)', line)
            if m:
                chromosomes[m.group(1)] = int(m.group(2))
                if logger:
                    logger.info('{}\t{}bp'.format(m.group(1), m.group(2)))
    proc.stdout.close()
    proc.wait()
    return chromosomes

def phi_correlation(c00, c01=None, c10=None, c11=None):
    if c01 is None:
        c00, c01, c10, c11 = c00[0:4]
    n0_ = c00 + c01
    n1_ = c10 + c11
    n_0 = c00 + c10
    n_1 = c01 + c11
    den = n0_ * n1_ * n_0 * n_1
    if den == 0:
        phi = 0
    else:
        try:
            phi = (c00 * c11 - c01 * c10) / np.sqrt(den)
        except:
            # print(c00, c01, c10, c11, n0_, n1_, n_0, n_1)
            raise
    return phi

def get_ordered_leaves(matrix, metrics=None):
    if matrix.shape[0] != matrix.shape[1]:
        return get_ordered_leaves_nondiagonal(matrix, metrics)
    import scipy.cluster.hierarchy
    # print(matrix)
    n = matrix.shape[0]
    X = []
    for i in range(n):
        for j in range(i + 1, n):
            X.append(matrix[i,j])
    Z = scipy.cluster.hierarchy.ward(X)
    ll = scipy.cluster.hierarchy.leaves_list(scipy.cluster.hierarchy.optimal_leaf_ordering(Z, matrix))
    return ll

def get_ordered_leaves_nondiagonal(matrix, metrics=None): # non-diagonal
    import scipy.cluster.hierarchy
    # print(matrix)
    n = matrix.shape[0]
    X = []
    dmat = np.zeros((n, n))
    if metrics is None:
        metrics = lambda x, y: np.dot(x-y, x-y)
    for i in range(n):
        vi = matrix[i,:]
        dmat[i,i] = 0
        for j in range(i + 1, n):
            vj = matrix[j,:]
            # print(vi, vj)
            dmat[i,j] = dmat[j,i] = metrics(vi, vj)
    return get_ordered_leaves(dmat)
    # Z = scipy.cluster.hierarchy.ward(X)
    # ll = scipy.cluster.hierarchy.leaves_list(scipy.cluster.hierarchy.optimal_leaf_ordering(Z, matrix))
    # return ll

