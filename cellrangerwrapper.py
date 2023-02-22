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
import tkutil
# def __get_logger(logger=None, logfile=None):
#     if logger is None:
#         logger = logging.getLogger(sys._getframe().f_code.co_name)
#     # set logging
#     def _set_log_handler(logger, handler):#, verbose):
#         handler.setFormatter(logging.Formatter('%(asctime)s %(name)s:%(lineno)s %(funcName)s [%(levelname)s]: %(message)s'))
#         logger.addHandler(handler)
#         return logger
#     _set_log_handler(logger, logging.StreamHandler())
#     if logfile is not None:
#         _set_log_handler(logger, logging.FileHandler(logfile))
#     # logger.setLevel(logging.ERROR)
#     logger.propagate = False
#     return logger

# def instanciate_standard_logger(name=None):
#     if name is None: name = __name__
#     logger = logging.getLogger(name)
#     return __get_logger(logger)

def _load_items(dirname, **kwargs):
    name = kwargs.get('name')
    column = kwargs.get('column', -1)
    trim_suffix = kwargs.get('trim', False)
    fbz = os.path.join(dirname, f'{name}.tsv.gz')
    fb = os.path.join(dirname, f'{name}.tsv')
    items = []
    if os.path.exists(fbz):
        with gzip.open(fbz) as fi:
            for line in fi:
                items.append(line.decode('utf-8').strip())
    else:
        with open(fb) as fi:
            for line in fi:
                items.append(line.strip())
    if column >= 0:
        data = []
        for line in items:
            data.append(line.split('\t')[column])
        items = data
    if trim_suffix:
        data = []
        for line in items:
            data.append(re.split('\\W', line)[0])
        items = data
    return items

def load_barcodes(dirname, **kwargs):
    """Load barcodes.tsv or barcodes.tsv.gz"""
    kwargs['name'] = 'barcodes'
    return _load_items(dirname, **kwargs)

def load_features(dirname, **kwargs):
    kwargs['name'] = 'features'
    return _load_items(dirname, **kwargs)

def load_sparse_matrix(dirname:str, **kwargs):
    """Load matrx.mtx
    """
    import gzip
    fm = os.path.join(dirname, 'matrix.mtx')
    mtz = os.path.join(dirname, 'matrix.mtx.gz')

    if os.path.exists(mtz):
        mtx = scipy.io.mmread(mtz)
    elif os.path.exists(fm):
        mtx = scipy.io.mmread(fm)
    else:
        raise Exception('{} does not include data'.format(dirname))

    return mtx

def count_total_reads(dirname:str, **kwargs):
    """bycell : """
    import hashlib, pickle
    logger = kwargs.get('logger', logging.getLogger())
    bycell = kwargs.get('bycell', False)
    bygene = kwargs.get('bygene', False)
    if bycell is False and bygene is False:
        bycell = True
    forced = kwargs.get('forced', False)
    if bycell == bygene:
        raise Exception('select bycell or bygene')
    md5 = hashlib.md5()
    md5.update(os.path.abspath(dirname).encode('ascii'))
    fm = os.path.join(dirname, 'matrix.mtx')
    mtz = os.path.join(dirname, 'matrix.mtx.gz')
    fn_cache = os.path.join(dirname, '.total.count.cache')
    
    if os.path.exists(mtz):
        dt = os.stat(mtz).st_mtime
        # md5.update('st_mtime={:.0f}'.format(os.stat(mtz).st_mtime).encode('ascii'))
        # flag = tkutil.check_file((mtz, dtst))
        fm = mtz
    elif os.path.exists(fm):
        dt = os.stat(mtz).st_mtime
        # md5.update('st_mtime={:.0f}'.format(os.stat(fm).st_mtime).encode('ascii'))
        # dt = os.stat(mtz).st_mtime).encode('ascii')
        
    logger.info('cache={}, matrix={}, dt={}'.format(fn_cache, fm, dt))
    if tkutil.check_file(fn_cache, dt) and not forced:
        logger.info('loading from cache')
        with gzip.open(fn_cache) as fi:
            cached = pickle.load(fi)
            if bycell:
                return cached['bycell']
            elif bygene:
                return cached['bygene']

    logger.info('load from {}'.format(fm))
    barcodes = load_barcodes(dirname)
    features = load_features(dirname)
    
    # fn_cache = os.path.join(dirname, '.{}.total.cache'.format(md5.hexdigest()[0:8]))
    # if os.path.exists(fn_cache) and os.path.getsize(fn_cache) > 100:
    #     logger.info('loading total counts from cache')
    #     return pd.read_csv(fn_cache, sep='\t', dtype=np.int32)
    logger.info('loading matrix to count total reads')
    mtx = load_sparse_matrix(dirname)    
    # bycell = np.array(mtx.sum(axis=0).tolist()[0]).reshape(-1, 1)
    # bygene= np.array(mtx.sum(axis=1).tolist()[0]).reshape(-1, 1)
    
    df_cell = pd.DataFrame(np.array(mtx.sum(axis=0)).reshape(-1, 1), index=barcodes, columns=['count', ])
    # df.to_csv(fn_cache, sep='\t')
    df_gene = pd.DataFrame(np.array(mtx.sum(axis=1)).reshape(-1, 1), index=features, columns=['count', ])
    data = {'bycell':df_cell, 'bygene':df_gene}
    # data['bycell'] = df_cell
    # data['bygene'] = df_gene
    logger.info('saving count matrices to {}'.format(fn_cache))
    with gzip.open(fn_cache, 'wb') as fo:
        pickle.dump(data, fo)
    if bycell:
        return df_cell
    else:
        return df_gene
    

def load_reads_from_sparse_matrix(srcdir:str, **kwargs)->pd.DataFrame:
    verbose = kwargs.get('verbose', False)
    fn_cache = os.path.join(srcdir, '.count.cache')
    if os.path.exists(fn_cache) and os.path.getsize(fn_cache) > 1000:
        df = pd.read_csv(fn_cache, sep='\t', dtype=np.int32)
        return df
    mtx = load_sparse_matrix(srcdir)
    s = np.sum(mtx, axis=0).tolist()[0]
    mtx[mtx>1] = 1
    t = np.sum(mtx, axis=0).tolist()[0]
    del mtx
    df = pd.DataFrame([s, t], columns=['n_Reads', 'n_Features'], index=barcodes).T
    df.to_csv(fn_cache, sep='\t')
    return df
    # return counts_per_cell, features_per_cell

def load_counts(srcdir:str, **kwargs)->pd.DataFrame:
    if os.path.isfile(srcdir):
        return pd.read_csv(srcdir, sep='\t', index_col=0)
    else:
        return load_reads_from_sparse_matrix(srcdir)

def load_tpm(filename:str, **kwargs):
    output = kwargs.get('output', None)
    use_log = kwargs.get('use_log', None)
    verbose = kwargs.get('verbose', False)
    forced = kwargs.get('forced', False)

    if output is None:
        if use_log:
            output = os.path.join(os.path.dirname(filename), '.{}.logtpm'.format(os.path.basename(filename)))
        else:
            output = os.path.join(os.path.dirname(filename), '.{}.tpm'.format(os.path.basename(filename)))
    if filename.endswith('.tpm') or (os.path.exists(output) and os.path.getsize(output) > 1000):
        if not forced:
            return output
    if verbose:
        sys.stderr.write('converting counts to TPM {}\n'.format(output))
    t = pd.read_csv(filename, index_col=0, sep='\t')
    tpm = t / np.array(np.sum(t, axis=0)) * 1e6
    if use_log:
        tpm = np.log2(tpm + 1)
    tpm.to_csv(output, sep='\t', float_format='%.2f')
    return output

def save_sparse_matrix(dstdir, matrix, barcodes=None, features=None, logger=None):
    import scipy.io
    import scipy.sparse
    import io
    import subprocess
    import pandas
    os.makedirs(dstdir, exist_ok=1)
    fn_barcode = os.path.join(dstdir, 'barcodes.tsv.gz')
    fn_mtx = os.path.join(dstdir, 'matrix.mtx')
    fn_features = os.path.join(dstdir, 'features.tsv.gz')
    if logger is None:
        logger = logging.getLogger(__name__)

    if isinstance(matrix, pandas.DataFrame):
        barcodes = matrix.columns
        features = matrix.index
        matrix = scipy.sparse.coo_matrix(matrix.values)
    elif isinstance(matrix, np.ndarray) or isinstance(matrix, list):
        matrix = scipy.sparse.csr_matrix(matrix)
    if barcodes is None or features is None:
        raise Exception('features and barcodes are required')
    if matrix.shape[0] != len(features) or matrix.shape[1] != len(barcodes):
        raise Exception('inconsistent matrix size {}x{}, n_features={}, n_barcodes={}'.format(matrix.shape[0], matrix.shape[1], len(features), len(barcodes)))
    
    logger.info('saving sparse matrix to {}'.format(fn_mtx))
    scipy.io.mmwrite(fn_mtx, matrix)
    fn_mtx_gz = fn_mtx + '.gz'
    if os.path.exists(fn_mtx_gz): os.unlink(fn_mtx_gz)
    cmd = 'pigz', '-p', '4', fn_mtx
    logger.info('compressing {}'.format(fn_mtx))
    proc = subprocess.Popen(cmd).wait()

    with io.TextIOWrapper(gzip.open(fn_barcode, 'wb'), encoding='utf-8') as f1,\
        io.TextIOWrapper(gzip.open(fn_features, 'wb'), encoding='utf-8') as f2:
        for c in barcodes:
            f1.write(c + '\n')
        for g in features:
            f2.write(g + '\n')

def convert_sparse_matrix_to_count(srcdir:str, filename:str=None, **kwargs):
    """Convert sparse matrix to tsv"""
    verbose = kwargs.get('verbose', False)
    forced = kwargs.get('forced', False)
    filename = kwargs.get('filename', None)
    feature_field = kwargs.get('field', None)

    if filename is None:
        filename = os.path.join(srcdir, 'count.tsv')
        if feature_field is not None:
            filename = os.path.join(srcdir, 'count.{}.tsv'.format(feature_field))

    if os.path.exists(filename) and os.path.getsize(filename) > 0 and not forced:
        return filename
    if verbose:
        sys.stderr.write('\033[Kloading barcodes\r')
    barcodes = load_barcodes(srcdir)
    if verbose:
        sys.stderr.write('\033[Kloading features\r')
    features = load_features(srcdir)
    if feature_field is not None:
        f2i = {}
        n_features = 0
        genes = []
        for i, feature in enumerate(features):
            items = feature.split('\t')
            if len(items) <= feature_field:
                feature_field = 0
            gene = items[feature_field]
            if gene not in f2i: f2i[gene] = []
            genes.append(gene)
            f2i[gene].append(i)
            n_features += 1
        if n_features != len(f2i):
            if verbose:
                sys.stderr.write('degeneration required\n')
        else:
            fetures = genes
            f2i = None
    else:
        f2i = None

    if verbose:
        sys.stderr.write('\033[Kloading sparse matrix\r')
    mtx = load_sparse_matrix(srcdir)
    if f2i is not None:
        if verbose:
            sys.stderr.write('\033[Kaggregating {} features into {} rows\r'.format(mtx.shape[0], len(f2i)))
        mtx = mtx.tocsr()
        rows = []
        ft_ = []
        for gene in sorted(f2i.keys()):
            idx = f2i[gene]
            ft_.append(gene)
            if len(idx) == 1:
                row = mtx[idx[0], :]
            else:
                row = scipy.sparse.csr_matrix(np.sum(mtx[idx, :], axis=0))
            rows.append(row)
        mtx = scipy.sparse.vstack(rows)
        features = ft_

    if verbose:
        sys.stderr.write('\033[Ksaving count matrix to {}\n'.format(filename))
    df = pd.DataFrame(mtx.toarray(), columns=barcodes, index=features, dtype=np.int32)
    df.to_csv(filename, sep='\t')
    return filename

def _check_sparse_matrix(srcdir):
    for f in 'matrix.mtx', 'barcodes.tsv', 'features.tsv':
        flag = False
        for z in ('', '.gz'):
            fn = os.path.join(srcdir, f + z)
            if os.path.exists(fn) and os.path.getsize(fn) > 10:
                flag = True
                break
        if not flag:
            return False
    return True

def load_gene_expression(filename:str, genes:list, **kwargs)->pd.DataFrame:
    """Load expression DataFrame from count/tpm file or sparse matrix.
    If the genes are cached, this function returns DataFrame first.
    """
    import hashlib, sqlite3, gzip, json
    verbose = kwargs.get('verbose', False)
    forced = kwargs.get('forced', False)
    logger = kwargs.get('logger', None)
    if logger is None:
        logger = logging.getLogger(__name__)
        if verbose:
            logger.setLevel(logging.DEBUG)
    
    # print(filename)
    # print(os.path.isdir(filename))
    # print(filename, os.path.isdir(filename), _check_sparse_matrix(filename))
    # if os.path.isdir(filename):
    #     logger.info('loading sparse matrix')
    #     if _check_sparse_matrix(filename):
    #         feature_field = kwargs.get('feature_field', 1)
    #         filename = convert_sparse_matrix_to_count(filename, field=feature_field, verbose=verbose, forced=forced)
    #     else:
    #         raise Exception("invalid directory")

    if os.path.isfile(filename):
        srcdir = os.path.dirname(filename)
    else:
        srcdir = filename

    genes = set(genes)
    if isinstance(genes, str):
        genes = [genes, ]
    filename_db = os.path.join(srcdir, '.' + os.path.basename(filename) + '.expr.db')
    
    logger.info('expression database file is {}'.format(filename_db))
    expression = {}
    try:
        with sqlite3.connect(filename_db) as cnx:
            pass
    except:
        sys.stderr.write('cannot open ' + filename_db)
    with sqlite3.connect(filename_db) as cnx:
        cur = cnx.cursor()
        cur.execute('create table if not exists expression(gene not null primary key, data blob)')
        cached_genes = []
        genes_loaded = []
        if not forced:
            cstr = ''
            if len(genes) < 250:
                for gene in genes:
                    if cstr != '':
                        cstr += ' or '
                    cstr += 'gene="{}"'.format(gene)
                sqlcom = 'select gene, data from expression where ' + cstr
            else:
                sqlcom = 'select gene, data from expression'
            if verbose:
                logger.info(sqlcom)
            cur.execute(sqlcom)
            for r in cur.fetchall():
                gene = r[0]
                if gene in genes:
                    if r[1] is not None and len(r[1]) > 0:
                        # logger.info(f'reading {gene} from cache')
                        values = json.loads(gzip.decompress(r[1]).decode('utf-8'))
                        if len(values) > 0:
                            expression[gene] = values
                    else:
                        genes_loaded.append(gene)
        n_genes = len(genes)
        existing = set()
        cur.execute('select gene from expression')
        for r in cur.fetchall():
            existing.add(r[0])

    # load columns        
    if os.path.isfile(filename):
        with open(filename) as fi:
            columns = [x_.strip('"\n') for x_ in fi.readline().split('\t')]
            if columns[0] != '' and (columns[0].lower() not in ('tracking_id', 'gene', 'transcript_id', 'gene_id')):
                columns = columns
            else:
                columns = columns[1:]
    else:
        columns = load_barcodes(srcdir)

    # all data cached
    if len(expression) == n_genes:
        m_ = []
        g_ = []
        for g in sorted(expression.keys()):
            g_.append(g)
            m_.append(expression[g])
        return pd.DataFrame(np.array(m_, dtype=np.float32), index=g_, columns=columns)
    genes_to_load = set([g for g in genes if g not in genes_loaded])
    genes_to_save = []
    logger.info('{}'.format(','.join(genes_to_load)))
    if len(genes_to_load) > 0:
        logger.info('loading {} genes from {} : {}'.format(len(genes_to_load), filename, ','.join(genes_to_load)))
        if os.path.isdir(filename): # sparse matrix
            logger.info('load sparce matrix from ' + filename)
            srcdir = filename
            columns = load_barcodes(srcdir)
            features = load_features(srcdir)
            idx = []
            ft_ = []
            for i, f in enumerate(features):
                f_ = f.split('\t')[0]
                if f_ in genes_to_load:
                    idx.append(i)
                    ft_.append(f_)
            mat = load_sparse_matrix(srcdir).tocsr()[idx,:].toarray()
            for i, f_ in enumerate(ft_):
                # logger.info(f_)
                # logger.info(i)
                # logger.info(mat[i,:])
                # logger.info(f_)
            # for i, j in enumerate(idx):
                expression[f_] = mat[i,:].tolist()
                genes_to_save.append(f_)
        else:                
            with open(filename) as fi:
                columns = [x_.strip('"\n') for x_ in fi.readline().split('\t')]
                if columns[0] != '' and (columns[0].lower() not in ('tracking_id', 'gene', 'transcript_id', 'gene_id')):
                    columns = columns
                else:
                    columns = columns[1:]

                for line in fi:
                    gene, row = line[0:-1].split('\t', 1)
                    gene = gene.strip('"')
                    if gene in genes:
                        values = [float(x_) for x_ in row.split('\t')]
                        expression[gene] = values
                        genes_to_save.append(gene)
        with sqlite3.connect(filename_db) as cnx:
            cur = cnx.cursor()
            for gene in genes_to_save:
                e_ = expression[gene]
                valueobj = gzip.compress(json.dumps(e_).encode('utf-8'))
                if gene in existing:
                    cur.execute('update expression set data=? where gene=?', (valueobj, gene))
                else:
                    cur.execute('insert into expression (gene, data) values(?, ?)', (gene, valueobj))
            for gene in genes_to_load:
                if gene not in genes_to_save:
                    # zobj = gzip.compress('[]'.encode('utf-8'))
                    logger.info(f'{gene} is set as NULL')
                    # logger.info(zobj)
                    cur.execute('insert into expression (gene, data) values(?, NULL)', (gene, ))
            cnx.commit()
    m_ = []
    g_ = []
    for g in sorted(expression.keys()):
        g_.append(g)
        m_.append(expression[g])
    if len(m_) == 0:
        return pd.DataFrame([], columns=columns)
    return pd.DataFrame(np.array(m_, dtype=np.float32), index=g_, columns=columns)

def convert_multi_to_gex(arguments=None):
    import subprocess
    import scipy.io
    if arguments is None:
        paresr = argparse.ArgumentParser()
        parser.add_argument('-i', nargs='+')
        parser.add_argument('-o')
        parser.add_argument('--verbose', action='store_true')
        args = parser.parse_args()
    else:
        args = arguments
    logger = tkutil.instanciate_standard_logger(sys._getframe().f_code.co_name)
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    srcdir = args.i[0]
    outdir = args.o
    os.makedirs(outdir, exist_ok=1)
    gex_index = []

    logger.info('copy barcodes')
    barcodes = load_barcodes(srcdir)
    fn_bco = os.path.join(outdir, 'barcodes.tsv.gz')
    with gzip.open(fn_bco, 'wb') as fo:
        for b in barcodes:
            fo.write((b + '\n').encode('utf-8'))
        # fo.write(barcodes.encode('ascii'))

    logger.info('detecting gene index')
    
    fn_ft = os.path.join(srcdir, 'features.tsv')
    fn_ftz = os.path.join(srcdir, 'features.tsv.gz')
    fn_fto = os.path.join(outdir, 'features.tsv.gz')

    if os.path.exists(fn_ft):
        ist = open(fn_ft)
    elif os.path.exists(fn_ftz):
        ist = io.TextIOWrapper(gzip.open(fn_ftz), encoding='utf-8')
    else:
        logger.error(f'no feature file in {srcdir}')
        raise Exception(f'no feature file in {srcdir}')
    index = 0
    gex_index = []
    with gzip.open(fn_fto, 'wb') as fo:
        for line in ist:
            items = line.strip().split('\t')
            if line[0] in '%#!': continue
            index += 1
            if len(items) < 3 or items[2] == 'Gene Expression':
                gex_index.append(index)
                fo.write(line.encode('utf-8'))
        ist.close()

    logger.info('editing matrix')
    fn_mto = os.path.join(outdir, 'matrix.mtx')
    for fn in 'matrix.mtx', 'matrix.mtx.gz':
        fn_mt = os.path.join(srcdir, fn)
        if os.path.exists(fn_mt):
            logger.info(f'loading matrix from f{fn_mt}')
            mm = scipy.io.mmread(fn_mt).tocsr()[gex_index]
            scipy.io.mmwrite(fn_mto, mm, field='integer')
            if os.path.exists(fn_mto + '.gz'):
                os.unlink(fn_mto + '.gz')
            logger.info('compressing sparse matrix')
            subprocess.Popen(['pigz', '-p', '4', fn_mto]).wait()
            break

def merge_sparse_matrices(args=None):
    if args is None:
        parser = argparse.ArgumentParser()
        parser.add_argument('-i', nargs='+')
        parser.add_argument('-o')
        parser.add_argument('--verbose', action='store_true')
        pass
    srcdirs = args.i
    dstdir = args.o
    logger = tkutil.instanciate_standard_logger(sys._getframe().f_code.co_name)
    if args.verbose:
        logger.setLevel(logging.DEBUG)
        logger.info('logger')
    
    os.makedirs(dstdir, exist_ok=1)
    common_features = None
    candidates = ['.', 'outs/filtered_feature_bc_matrix', 'outs/raw_feature_bc_matrix', 'filtered_feature_bc_matrix', 'raw_feature_bc_matrix', ]
    merged_barcodes = []
    matrices = []
    titles = []
    logger.info('merge matrices')
    for i, srcdir in enumerate(srcdirs):
        title_ = title = os.path.basename(srcdir)
        count = 1
        while title_ in titles:
            title_ = f'{title}.{count}'
            count += 1
        titles.append(title_)
        title = title_
        basedir = None
        for candidate in candidates:
            dn = os.path.join(srcdir, candidate)
            if os.path.exists(dn) and os.path.isdir(dn) and os.path.exists(os.path.join(dn, 'matrix.mtx.gz')):
                basedir = dn
                break
        if basedir is None:
            logger.warn(f'no matrix in {srcdir}')
            continue
        barcodes = load_barcodes(basedir)
        features = load_features(basedir, column=-1)
        if common_features is None:
            common_features = features
        else:
            if len(common_features) != len(features):
                logger.error(f'size of features imcompatible in {srcdir}')
                continue
            compatible = True
            for i, f in enumerate(features):
                if f != common_features[i]:
                    logger.error(f'incompatible features in {srcdir}')
                    compatible = False
                    break
            if not compatible:
                continue
        logger.info('{} loading {}x{} matrix from {}'.format (title, len(features), len(barcodes), basedir))
        merged_barcodes += ['{}:{}'.format(title, b) for b in barcodes]
        mat = load_sparse_matrix(basedir)
        matrices.append(mat)
    logger.info('saving {}x{} matrix to {}'.format(len(features), len(merged_barcodes), dstdir))
    save_sparse_matrix(dstdir, scipy.sparse.hstack(matrices), merged_barcodes, common_features)
    
def edit_hdf(arguments=None):
    """Read molecule info
    """
    import h5py
    if arguments is not None:
        args = arguments
    else:
        parser = argparse.ArgumentParser()
        parser.add_argument('-i', metavar='')
        parser.add_argument('-o', default='h5')
        parser.add_argument('--barcode', metavar='tsv file')
        # parser.add_argument('--feature', metavar='tsv file')
        parser.add_argument('--verbose', action='store_true')
        args, cmds = parser.parse_known_args()  
    outdir = args.o
    verbose = args.verbose
    os.makedirs(outdir, exist_ok=1)
    fn_bcd = args.barcode
    logger = tkutil.instanciate_standard_logger(sys._getframe().f_code.co_name)
    if verbose: logger.setLevel(10)

    barcodes = []
    with open(fn_bcd) as fi:
        for line in fi:
            barcodes.append(re.sub('\\W\\d$', '', line.strip().split('\t')[0]))
    accepted_barcodes = set(barcodes)
    logger.info('{} barcodes accepted'.format(len(barcodes)))

    file_out = os.path.join(outdir, 'molecule_info.h5')
    # copy_fields = ['metrics_json', '', 'gem_group', 'library_info', ]

    features = ['features', 'feature_idx', ]
    barcodes = ['barcode_idx', 'barcode_info', 'barcodes']
    count = ['count']

    dataset_fields = [
        'barcode_idx', 
        'barcode_info/genomes',
        'barcode_info/pass_filter',
        'barcodes',
        'count',
        'feature_idx',
        'features/_all_tag_keys',
        'features/feature_type',
        'features/genome',
        'features/id',
        'features/name',
        'gem_group',
        'library_idx',
        'library_info',
        'metrics_json',
        'umi',
        'umi_type'
    ]
    file_in = args.i

    with h5py.File(file_in, 'r') as src, h5py.File(file_out, 'w') as dst:
        passed = src['/barcode_info/pass_filter']#barcodes/pass_filter']
        # print(passed.dtype)
        # print(passed.__class__)
        barcodes = src['barcodes']
        passfilter = []
        nextstep = 100
        for i, row in enumerate(passed):
            j = row[0]
            # print(i, j)
            b = re.sub('\\W\\d$', '', barcodes[j].decode('ascii'))
            if i > nextstep: 
                sys.stderr.write('\033[K {}\t{}\t{}\r'.format(i + 1, len(passfilter), b))
                nextstep += 100
                # break
            if b in accepted_barcodes:
                passfilter.append(row)
        logger.info('{} cells were selected from {} passed/{} total'.format(len(passfilter), len(passed), len(barcodes)))
        dst.create_dataset('barcode_info/pass_filter', data=np.array(passfilter, dtype=np.uint64))
        for attr in attributes:
            dst.create_attribute
        for field in dataset_fields:#src.keys():#copy_fields:
            logger.info(f'cloning {field}')
            if field =='barcode_info/pass_filter': 
                continue
            obj = src[field]
            # print(field, obj.__class__)
            dst.create_dataset(field, data=src[field])
        dst.flush()

def main():
    logger = logging.getLogger(sys._getframe().f_code.co_name)
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', nargs='+')
    parser.add_argument('-o')
    parser.add_argument('--verbose', action='store_true')

    args, carg = parser.parse_known_args()

    verbose = args.verbose
    if verbose:
        logger.setLevel(logging.DEBUG)
    cmd = carg[0]
    if cmd == 'multi2gex':
        convert_multi_to_gex()#args)
    elif cmd == 'merge':
        merge_sparse_matrices()
    elif cmd == 'edithdf':
        edit_hdf()#args)
    else:
        raise Exception(f'command {cmd} not supported')

if __name__ == '__main__':
    main()
