import os, sys, re
import pandas as pd
import numpy as np
import scipy.stats
import argparse
import statsmodels.stats.multitest
import pymultiscale.anscombe
import logging, argparse


"""Detection marker genes using ANOVA 

-i [filename] normalized tsv data of expressions
-c [filename] cluster.tsv
-o [filename] anova.tsv
--marker gene1 gene2 ... genes showing in process
--log log value
--cluster cluster1 cluster2 ...
--min-size minimum cluster size to consider
--add-ratio append nonzero percent columns
--seurat-dir output directory of Seurat analysis
--dispersion 1
--normalize-barcode skip barcode validation
"""

def __get_logger(logger=None, logfile=None):
    if logger is None:
        logger = logging.getLogger(sys._getframe().f_code.co_name)
    # set logging
    def _set_log_handler(logger, handler):#, verbose):
        handler.setFormatter(logging.Formatter('%(asctime)s %(name)s:%(lineno)s %(funcName)s [%(levelname)s]: %(message)s'))
        logger.addHandler(handler)
        return logger
    _set_log_handler(logger, logging.StreamHandler())
    if logfile is not None:
        _set_log_handler(logger, logging.FileHandler(logfile))
    # logger.setLevel(logging.ERROR)
    logger.propagate = False
    return logger

def save_loadings(matrix, barcodes=None, features=None, filename='loadings.tsv', n_components=2, n_genes=100, logger=None):
    if logger is None:
        logger = logging.getLogger()
    logger.info('calculating loadings')
    import sklearn.decomposition
    if isinstance(matrix, pd.DataFrame):
        barcodes = matrix.columns
        features = matrix.index
    else:
        if barcodes is None:
            barcodes = [str(i+1) for i in range(matrix.shape[1])]
        if features is None:
            features = [str(i+1) for i in range(matris.shape[0])]
    pca = sklearn.decomposition.IncrementalPCA(n_components=n_components)
    res = pca.fit_transform(matrix)
    loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
    table = []
    columns = []
    for i in range(n_components):
        columns += [f'PCA_{i+1}', f'Loading_{i+1}']
        gnames = []
        lvalues = []
        for j in sorted(range(loadings.shape[0]), key=lambda j_:loadings[j_,i] ** 2, reverse=True):
            gnames.append(features[j])
            lvalues.append(loadings[j,i])
        table.append(gnames)
        table.append(lvalues)
        if len(table) >= n_genes * 2:
            break
    df = pd.DataFrame(table, index=columns).T
    if filename.endswith('.xlsx'):
        df.to_excel(filename)
    else:
        df.to_csv(filename, sep='\t', float_format='%.5f')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', default=None, metavar='filename', help='cluster data')
    parser.add_argument('-e', default=None, metavar='filename', help='expression data')
    parser.add_argument('-o', '--outdir', help='output direcotry, default:out')
    # parser.add_argument('--cluster', nargs='+', default=None, metavar='clusters', help='use selected clusters')
    parser.add_argument('--marker', nargs='+', default=[])
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--normalize-barcode', action='store_true')
    parser.add_argument('--log', action='store_true')
    parser.add_argument('--use-vst', action='store_true')
    parser.add_argument('--min-size', type=int, default=100)
    parser.add_argument('--seurat-dir', default=None)
    parser.add_argument('--dispersion', default=0.1, type=float)
    parser.add_argument('--loadings', type=int, default=0, help='number of high loading genes, default 0 (not write)')
    # parser.add_argument('--pca', action='store_true', help='Dump PCA loadings using given conditions')
    # parser.add_argument('--add-ratio', action='store_true', help='Add percent of cells having signal')
    args = parser.parse_args()
    seurat_dir = args.seurat_dir
    outdir = args.outdir
    dispersion = args.dispersion
    ignore_barcode_validation = not args.normalize_barcode
    if seurat_dir is not None:
        fn_expr = os.path.join(seurat_dir, 'seurat.normalized.tsv')
        fn_clus = os.path.join(seurat_dir, 'seurat.clusters.tsv')
    else:
        fn_expr = fn_clus = None
    use_log = args.log
    verbose = args.verbose
    min_cluster_size = args.min_size
    if args.e is not None:
        fn_expr = args.e
    if args.c is not None:
        fn_clus = args.c
    if outdir.endswith('.tsv') or outdir.endswith('.csv') or outdir.endswith('.txt'):
        outdir = os.path.dirname(outdir)
    os.makedirs(outdir, exist_ok=1)
    fn_anova = os.path.join(outdir, 'anova.tsv')
    fn_loading = os.path.join(outdir, 'loading.tsv')
    use_vst = args.use_vst
    n_loadings = args.loadings
    cluster_titles = {}
    logger = __get_logger()
    if verbose:
        logger.setLevel(logging.DEBUG)

    if fn_clus is not None:
        clust = pd.read_csv(fn_clus, sep='\t', index_col=0)
        clusters = clust.iloc[:,0].values
        barcodes = [x_.strip('"').replace('-', '.') for x_ in clust.index]
    else:
        clus = {}
        clusters = []
        barcodes = []
        with open(fn_expr) as fi:
            for item in fi.readline().split('\t'):
                item = item.strip('"\n"')
                if item == '': continue
                m = re.match('(.*?)\\D\\d+$', item)                    
                if m:
                    barcodes.append(m.group(0))
                    cluster_name = m.group(1)
                    if cluster_name not in clus:
                        cluster_titles[len(clus)] = cluster_name
                        clus[cluster_name] = len(clus)
                    clusters.append(clus[cluster_name])
                else:
                    clusters.append(-1)
                    barcodes.append(item)
        # print(cluster_titles)

    n_clusters = max(clusters) + 1
    group2index = [[] for i in range(n_clusters)]
    barcode2index = {}
    index2barcode = {}
    barcode2cluster = {}
    unclassified = set()
    markerset = set(args.marker)

    for i, c in enumerate(clusters):
        barcode = barcodes[i]
        if ignore_barcode_validation:
            barcode = barcodes[i]
        else:
            barcode = re.sub('[^ACGT]', '', barcodes[i])
            # barcode = re.sub('^X(\\d+.*$)', '\\1', barcodes[i])
            # barcode = re.sub('[:\\-]', '.', barcode)
        if c >= 0:
            barcode2cluster[barcode] = c
            barcode2index[barcode] = i
        else:
            unclassified.add(barcode)

    def get_available_clusters(group2index, min_cluster_size):
        available_clusters = []
        if min_cluster_size > 0:
            for c in range(n_clusters):
                print(c, len(group2index[c]))
                if len(group2index[c]) >= min_cluster_size:
                    available_clusters.append(c)
        else:
            available_clusters = list(range(n_clusters))
        return available_clusters

    # print(unclassified)
        # if verbose:
        #     sys.stderr.write('barcodes={}, cluster={}, unclassified={}\n'
        #     .format(len(barcode2index), len(available_clusters), len(unclassified)))
    results = []
    # signals = []
    pvalues = []

    if os.path.isdir(fn_expr):
        srcdir = fn_expr
        import cellrangerwrapper    
        logger.info('loading sparse matrix from {}'.format(srcdir))
        barcodes = cellrangerwrapper.load_barcodes(srcdir)
        # print(len(barcodes))
        # print(barcodes[0:10])
        # print(list(barcode2cluster.keys())[0:10])
        # # barcodes = [re.sub('^X(\\d+.*)', '\\1', x_) for x_ in barcodes]
        features = cellrangerwrapper.load_features(srcdir)
        for i, item in enumerate(barcodes):
            if not ignore_barcode_validation:
                barcode = re.sub('^X(\\d+.*$)', '\\1', barcodes[i])
                barcode = re.sub('[:\\-]', '.', barcode)
            # print(barcode)
            if barcode not in barcode2cluster:
                print(f'{barcode} is absent')
                continue
            barcode2index[barcode] = i
            # print(barcode, i)
            group2index[barcode2cluster[barcode]].append(i)
        available_clusters = get_available_clusters(group2index, min_cluster_size)
        print(available_clusters)
        exit()
        matrix = cellrangerwrapper.load_sparse_matrix(srcdir).tocsr()

        if n_loadings > 0:
            logger.info('calculating PCA loading')
            if use_vst:
                lmat = pymultiscale.anscombe.anscombe(matrix)
            else:
                lmat = matrix
            save_loadings(lmat, barcodes, features, filename=fn_loading, n_components=n_pca_components, logger=logger, )
        for i in range(matrix.shape[0]):
            row = matrix[i,:].toarray().reshape(-1)
            if use_log:
                values = np.log2(row + dispersion)
            elif use_vst:
                values = pymultiscale.anscombe.anscombe(row)
            else:
                values = row
            valset = []
            for cn in available_clusters:#group2index:
                try:
                    idx =np.array(group2index[cn])
                    valset.append(values[idx])
                except:
                    sys.stderr.write(str(idx) + '\n')
                    sys.stderr.write(str(values) + '\n')
                    raise
            means = [np.mean(vals) for vals in valset]
            if np.max(means) > 0:
                pvalue = scipy.stats.f_oneway(*valset).pvalue
            else:
                pvalue = 1.0
            gene = features[i]
            if pvalue < 0.01:
                ostr = gene
                for i, m in enumerate(means):
                    ostr += '\t{:.3f}'.format(m)
                sys.stderr.write('\033[K{}\r'.format(ostr))
                if gene in markerset:
                    sys.stderr.write('\n')
            pvalues.append(pvalue)
            results.append((gene, means, pvalue))
    else:        
        logger.info('loading tsv file {}'.format(fn_expr))
        if n_loadings > 0:
            logger.info('loading csv to calculate loading')
            matrix = pd.read_csv(fn_expr, sep='\t')
            if use_vst:
                lmat = pymultiscale.anscombe.anscombe(matrix)
            else:
                lmat = matrix
            save_loadings(matrix=lmat, filename=fn_loading, n_components=n_loadings, logger=logger)
            table = None
        # table = None if calculate_loadings is False else []
        table = None
        # print(barcode2index)
        with open(fn_expr) as fi:
            header = [re.sub('[^ACGT]', '', x_) for x_ in fi.readline().strip().split('\t')]
            if ignore_barcode_validation:
                for i in range(len(barcode2cluster)):
                    # print(i)
                    barcode2index[i] = i
                    group2index[barcode2cluster[i]].append(i)
                available_clusters = get_available_clusters(group2index, min_cluster_size)
            else:
                for i, item in enumerate(header):
                    barcode = re.match('[\\w+\\.\\-]+', item.strip('"')).group(0)
                    # barcode = barcode.replace('.', '-')
                    barcode = re.sub('^X(\\d+.*)', '\\1', barcode)
                    if barcode not in barcode2cluster:#index:
                        sys.stderr.write('{} is absent in cluster\n'.format(barcode))
                        continue
                    barcode2index[barcode] = i
                    c_ = barcode2cluster[barcode]
                    group2index[c_].append(i)
                    available_clusters = get_available_clusters(group2index, min_cluster_size)
                while 1:
                    line = fi.readline()
                    if line == '': break
                    items = line.strip().split('\t')
                    gene = items[0].strip('"')
                    values = np.array([float(x_) for x_ in items[1:]])
                    if use_log:
                        values = np.log2(values + dispersion)
                    elif use_vst:
                        values = pymultiscale.anscombe.anscombe(values)
                    valset = []
                    # signal = []
                    for cn in available_clusters:#group2index:
                        idx =np.array(group2index[cn])
                        valset.append(values[idx])
                    means = [np.mean(vals) for vals in valset]
                    if np.max(means) == np.min(means):
                        pvalue = 1.0
                    else:
                        pvalue = scipy.stats.f_oneway(*valset).pvalue
                    if pvalue < 0.01:
                        ostr = gene
                        for i, m in enumerate(means):
                            ostr += '\t{:.3f}'.format(m)
                        sys.stderr.write('\033[K{}\r'.format(ostr))
                        if gene in markerset:
                            sys.stderr.write('\n')
                    pvalues.append(pvalue)
                    results.append((gene, means, pvalue))

        if verbose:
            sys.stderr.write('\033[K\r')

    import statsmodels.stats.multitest
    mres = statsmodels.stats.multitest.multipletests(pvalues)
    n_genes = len(pvalues)
    # add_ratio = args.add_ratio        
    logger.info('saving results to {}'.format(fn_anova))
    with open(fn_anova, 'w') as fo:
        header = 'Gene\tMax\tMin\tpVal\tqVal'
        # header_ratio = ''
        for cn in available_clusters:
            header += '\t{}_{}'.format(cluster_titles.get(cn, 'C{}'.format(cn+1)), len(group2index[cn]))
            # header_ratio += '\t{}_pos'.format(cluster_titles.get(cn, 'C{}'.format(cn+1)))
        # if add_ratio:
        fo.write(header + '\n')
        # else:
        #     fo.write(header + '\t' + header_ratio + '\n')

        for i in sorted(range(n_genes), key=lambda i_:results[i_][0]):
            r = results[i]
            gene, means, pvalue = r
            minidx = np.argmin(means)
            maxidx = np.argmax(means)
            maxval = np.max(means)
            minval = np.min(means)
            if maxval == minval:
                minclus = maxclus = '.'
                pvalue = qvalue = 1.0
                # logfc = 0.0
            else:
                minclus = 'C{}'.format(available_clusters[minidx] + 1)
                maxclus = 'C{}'.format(available_clusters[maxidx] + 1)
                qvalue = mres[1][i]
                # logfc = (maxval - minval) if use_log else np.log2(maxval / minval)
            ostr = '{}\t{}\t{}\t{:.3f}\t{:.3e}'.format(gene, maxclus, minclus, pvalue, qvalue)
            for m in means:
                ostr += '\t{:.3f}'.format(m)
            # if add_ratio:
            #     signal = signals[i]
            #     for s in signal:
            #         ostr += '\t{:.3f}'.format(s)
            fo.write(ostr + '\n')


if __name__ == '__main__':
    main()