"""
Visualization script used in Ohinata et al. 2021.
To detect clusters of dead cells or those of multiplet cells, 
'count' command was applied and read count distribution was estimated.
Contaminated feeder cells weere removed by gene expression specific to 
them such as Ctcf using 'marker' command.
THe command 'cluster' displays colored clusters in 2d plot.
"""

import os, sys, re
import pandas as pd
import numpy as np
import scipy.io
import argparse
import umap
import plotly
import plotly.graph_objs as go
import plotly.subplots
import scipy.sparse
import pickle
import gzip
import logging
import tkutil
from cellrangerwrapper import *
import logging

def calculate_coordinate_color(pos, cmin, cmax, random_factor=0.4, colormode=0):
    deg = [(pos[i] - cmin[i])/(cmax[i] - cmin[i]) for i in range(len(pos))]
    # print(pos, cmin, cmax, deg)
    if len(deg) == 2:
        if colormode == 0:
            red = (1-deg[1])*(1.5-deg[0])
            blue = deg[0]
            green = deg[1] * .8
        elif colormode == 1:
            red = 1 - deg[0]
            green = 1 - deg[1]
            blue = deg[1]
        elif colormode == 2:
            red = deg[0]
            blue = (1-deg[0]) * (1-deg[1])
            green = deg[1]
        elif colormode == 3:
            red = (0.5-deg[0]) * 2 + deg[1] * .05
            green = (deg[0] - .5) * 2 - deg[1] * .1
            blue = 1 - abs(deg[0] - 0.5) * 4 - deg[1] * .15
        elif colormode == 4:
            red = (1-deg[0]) * (1-deg[1])
            blue = deg[0] * (1-deg[1])
            green = deg[1] - deg[0] * .1
        elif colormode == 5:
            red = deg[0] * deg[1]
            blue = 1 - deg[1] + deg[0] * .05
            green = (1-deg[0]) ** 2 + deg[1] * .1
        elif colormode == 6:
            green = deg[0] + (deg[1] - .6) * 2.5
            red = deg[0] * (1-deg[1])
            blue = (1-deg[0]) ** 2 - deg[1] * .1
        else:
            return None        
    else:
        red = deg[0]
        blue = deg[1]
        green = deg[2]
    values = np.array((red, green, blue))
    if random_factor > 0:
        values += (np.random.rand(3) - 0.5) * random_factor 
    rgb = np.array(values * 256, dtype=np.int32)
    rgb[rgb<0] = 0
    rgb[rgb>255] = 255
    return 'rgb({},{},{})'.format(rgb[0], rgb[1], rgb[2])

def load_reads_from_sparse_matrix(srcdir:str, **kwargs):#->pd.DataFrame:
    """
    Load total reads, MT reads and number of features from sparse matrix
    """
    verbose = kwargs.get('verbose', False)
    fn_cache = os.path.join(srcdir, '.count.cache')
    columns = ['n_Reads', 'n_Features', 'n_MT']
    if os.path.exists(fn_cache) and os.path.getsize(fn_cache) > 1000:
        df = pd.read_csv(fn_cache, sep='\t', index_col=0)
        flag_compatible = False
        if len(df.columns) == len(columns):
            flag_compatible = True
            for i, c in enumerate(df.columns):
                if c != columns[i]:
                    flag_compatible = False
                    break
        if flag_compatible:
            return df
    mtx = load_sparse_matrix(srcdir).tocsr()
    s = np.array(np.sum(mtx, axis=0).reshape(-1, 1))
    barcodes = load_barcodes(srcdir)
    ft = load_features(srcdir)
    midx = []
    for i, f in enumerate(ft):
        if f.lower().startswith('mt-'):
            midx.append(i)
    total_reads = np.asarray(mtx.sum(axis=0)).reshape(-1)
    mt_reads = np.asarray(mtx[midx,:].sum(axis=0)).reshape(-1)
    mtx[mtx>0] = 1
    n_features = np.asarray(mtx.sum(axis=0)).reshape(-1)
    df = pd.DataFrame(np.concatenate([s, total_reads, mt_reads], axis=1), columns=columns, index=barcodes)
    df.to_csv(fn_cache, sep='\t')
    return df
    # return counts_per_cell, features_per_cell
def load_reads_from_count_matrix(filename):
    # verbose = kwargs.get('verbose', False)
    srcdir = os.path.dirname(filename)
    print(srcdir)
    fn_cache = os.path.join(srcdir, '.count.from_tsv.cache')
    print(fn_cache)
    if os.path.exists(fn_cache) and os.path.getsize(fn_cache) > 1000:
        # print(fn_cache)
        df = pd.read_csv(fn_cache, sep='\t', index_col=0)
        return df
    with open(filename) as fi:
        barcodes = fi.readline().strip().split('\t')
        n_cells = len(barcodes)
        tbl = np.zeros((n_cells, 2), dtype=np.int32)
        for line in fi:
            items = line.strip().split('\t')
            for i, item in enumerate(items[1:]):
                if item != '0':
                    c = int(item)
                    tbl[i,0] += c
                    tbl[i,1] += 1
        pass
    df = pd.DataFrame(tbl, columns=['n_Reads', 'n_Features'], index=barcodes)
    df.to_csv(fn_cache, sep='\t')
    return df


def _load_sparse_matrix(dirname:str, **kwargs):
    """Load matrx.mtx
    """
    import gzip
    # fgz = os.path.join(dirname, 'features.tsv.gz')
    fm = os.path.join(dirname, 'matrix.mtx')
    mtz = os.path.join(dirname, 'matrix.mtx.gz')

    if os.path.exists(fm):
        mtx = scipy.io.mmread(fm)
    elif os.path.exists(mtz):
        mtx = scipy.io.mmread(mtz)
    else:
        raise Exception('{} does not include data'.format(dirname))

    return mtx

def draw_swarmplot():
    import collections
    import seaborn as sns
    import matplotlib.pyplot as plt
    import matplotlib
    import matplotlib.backends.backend_pdf
    # parser = argparse.ArgumentParser()
    # parser.add_argument('--verbose', action='store_true')
    # parser.add_argument('--seurat-dir')
    # parser.add_argument('-e')
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', help='expression tsv or sparse matrix')
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('-g', nargs='+')
    parser.add_argument('-c', help='cluster tsv file')
    parser.add_argument('-o', default='swarmplot')
    parser.add_argument('--seurat-dir', default='seurat', metavar='directory')
    parser.add_argument('--percentile', default=0.8, type=float, help='expression threshold, 0.9 as 90% in percentile mode')
    parser.add_argument('--open', action='store_true')
    parser.add_argument('--color')
    parser.add_argument('--ncol', type=int, default=3)
    parser.add_argument('--mode', choices=['swarm', 'strip'], default='swarm')
    parser.add_argument('--cluster', nargs='+')
    args, cmds = parser.parse_known_args()
    
    fn_pdf = args.o

    logger = tkutil.get_logger(os.path.basename(__file__))    
    if args.verbose:
        logger.setLevel(10)

    mode = args.mode
    auto_open = args.open
    basedir = args.seurat_dir
    logger.info('loading from ' + basedir)
    fn_expr = args.e or os.path.join(basedir, 'seurat.normalized.tsv')
    fn_clu = args.e or os.path.join(basedir, 'seurat.clusters.tsv')
    logger.info(basedir)
    logger.info(fn_expr)
    logger.info(fn_clu)

    genes = args.g
    ncol = args.ncol

    logger.info('loading expression from {}'.format(fn_expr))
    expr = load_gene_expression(fn_expr, genes=genes, logger=logger)
    # print(expr)
    logger.info('{} genes requested, {} genes loaded'.format(len(genes), expr.shape[0]))
    genes = [g for g in genes if g in expr.index]
    n_genes = expr.shape[0]
    
    if ncol <= 0:
        ncol = np.ceil(np.sqrt(n_genes))
    nrow = (n_genes + ncol - 1) // ncol
    logger.info('chart size col={}, row={} for {} genes'.format(ncol, nrow, n_genes))

    # col = row = 1
    # for i, gene in enumerate(genes):
    #     col += 1
    #     if col > ncol:
    #         row += 1
    #         col = 1
    
    logger.info('setup {} x {} graphhs'.format(nrow, ncol))
    sns.axes_style('white')
    sns.set(rc={'figure.figsize':(ncol * 4 + 1, nrow * 4.5+ 1), 'axes.facecolor':'white', 'axes.edgecolor':'black', 'axes.facecolor':'white'})
    doc = matplotlib.backends.backend_pdf.PdfPages(fn_pdf)
    fig, axes = plt.subplots(nrow, ncol)

    # clusters    
    df_clu = pd.read_csv(fn_clu, sep='\t', index_col=0)#DataFrame(columns=clusters)
    df_clu.columns = ['Cluster', ]
    n_clusters = np.max(df_clu.values) + 1
    n_cells = df_clu.shape[0]
    clusters = collections.OrderedDict()
    corder = np.full(n_cells, -1, dtype=np.int32)

    if args.cluster is None:
        for i in range(n_clusters):
            clusters[i] = []
    else:
        for c in args.cluster:
            cnum = -1
            if c.isdigit():
                cnum = int(c)
            else:
                m = re.match('C(\\d+)', c)
                if m:
                    cnum = int(m.group(1)) - 1
            if cnum < 0:
                logger.warninig('cannot resolve cluster name {}'.format(c))
                continue
            clusters[cnum] = []
    c2n = {}
    for c in clusters.keys():
        c2n[c] = len(c2n)
    idx_clustered = []
    cmat = []
    # print(clusters)
    for i, cn in enumerate(df_clu.values.reshape(-1)):
        if cn in clusters:
            clusters[cn].append(i)
            idx_clustered.append(i)
            cmat.append(['C{}'.format(cn + 1), c2n[cn]])
        else:
            cmat.append(['unc', -1])
    graph_data = []
    data = pd.concat([expr.T, pd.DataFrame(cmat, columns=['Cluster', '_cluster_order_'], index=expr.columns), ], axis=1)
    data = data.sort_values(by='_cluster_order_')
    data = data[data['_cluster_order_']>=0]
    
    import warnings
    warnings.simplefilter('ignore', category=UserWarning)
    
    colors = ['navy', ]
    col = row = 0    
    # print(axes)
    for gene in genes:
        logger.info('{}, {}\t{}'.format(row, col, gene))
        if mode == 'swarm':
            if data.shape[0] > 2000:
                try:
                    sampled = data.sample(2000).sort_index()
                except:
                    pass
            else:
                sampled = data
            ax_ = sns.swarmplot(data=sampled, x='Cluster', y=gene, palette=colors, size=2,
                edgecolor='black', marker='o', ax=axes[row][col], zorder=1).set_title(gene)
        else:
            ax_ = sns.stripplot(data=data, x='Cluster', y=gene, palette=colors, size=1.5, jitter=0.25, dodge=True,
                edgecolor='black', marker='o', ax=axes[row][col], zorder=1).set_title(gene)
        col += 1
        if col >= ncol:
            col = 0
            row += 1
    fig.tight_layout(pad=5.0)
    doc.savefig(pad_inches=1, bbox_inches='tight')
    doc.close()
    plt.clf()
    try:
        if auto_open:
            subprocess.Popen(['open', fn_pdf]).wait()
    except:
        pass
    return fn_pdf    
    
    # parser.add_argument('-s', '--sep', default='_', help='delimiter between name and barcode')
    # parser.add_argument('-u', help='coordinates such as UMAP or tSNE')
    # parser.add_argument('-o', default='chart')
    # parser.add_argument('--open', action='store_true')
    # args, cmds = parser.parse_known_args(arguments)
    
    
    
    
    
    """
    return filename
    """
    
def display_cell_distribution(arguments=None):
    """
    return {
        'scatter.html':filename, 
    }
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('-s', '--sep', default='_', help='delimiter between name and barcode')
    parser.add_argument('-u', help='coordinates such as UMAP or tSNE')
    parser.add_argument('-o', default='chart')
    parser.add_argument('--open', action='store_true')
    args = parser.parse_known_args(arguments)[0]

    logger = tkutil.instanciate_standard_logger(sys._getframe().f_code.co_name)
    verbose = args.verbose
    outdir = args.o
    os.makedirs(outdir, exist_ok=True)
    fn_chart = os.path.join(outdir, 'barcodes.html')
    if verbose: logger.setLevel(logging.DEBUG)
    auto_open = args.open
    fn_coord = args.u
    coord = pd.read_csv(fn_coord, sep='\t', index_col=0)
    barcodes= coord.index
    group2barcode = {}
    sep = args.sep
    for i, barcode in enumerate(coord.index):
        items = barcode.split(sep)
        gn = items[0]
        if gn not in group2barcode:
            group2barcode[gn] = []
        group2barcode[gn].append(i)
    logger.info(group2barcode.keys())
    max_clusters = 32
    others = []
    if len(group2barcode) > max_clusters:
        nums = sorted([len(bc) for bc in group2barcode.values()], reverse=True)
        thr_count = nums[max_clusters]
        for gn, idx in group2barcode.items():
            if len(idx) <= thr_count:
                others += idx
            group2barcode[gn] = []
    traces = []
    if len(others) > 0:
        xy = coord.values[others, 0:2]
        traces.append(go.Scattergl(x=xy[:,0], y=xy[:,1], name='Others', mode='markers', text=barcodes[others],
        marker=dict(size=2, color='gray')))
    for gn in sorted(group2barcode.keys()):
        idx = group2barcode[gn]
        if len(idx) > 0:
            xy = coord.values[idx, 0:2]
            traces.append(go.Scattergl(x=xy[:,0], y=xy[:,1], name=gn, mode='markers', text=barcodes[idx],
        marker=dict(size=4)))
    fig = go.Figure(traces)
    fig.update_layout(plot_bgcolor='white')
    fig.update_xaxes(linewidth=2, linecolor='black', showgrid=False, showline=True, mirror='ticks')
    fig.update_yaxes(linewidth=2, linecolor='black', showgrid=False, showline=True, mirror='ticks')
    plotly.offline.plot(fig, filename=fn_chart, auto_open=auto_open)
    return {
        'scatter.html':fn_chart,
    }

def show_marker(arguments=None):
    """
    return {'scatter.html':fn_scatter,
        'info':fn_info,
        'swarmplot.html':fn_swarm,
        'violin.html':fn_violin,
        'total.html':fn_total   
    }
    """
    if arguments is None:
        parser = argparse.ArgumentParser()
        parser.add_argument('-e', help='expression tsv or sparse matrix')
        parser.add_argument('--verbose', action='store_true')
        parser.add_argument('-g', nargs='+')
        parser.add_argument('-c', help='cluster')
        parser.add_argument('-u', help='coordinates such as UMAP or tSNE')
        parser.add_argument('--dim', type=int, default=2, choices=[2, 3], help='2D or 3D chart')
        parser.add_argument('-o', default='exprchart')
        parser.add_argument('--seurat-dir', default='seurat', metavar='directory')
        # parser.add_argument('--threshold', default=1, type=float, help='expression threshold, 0.9 as 90% in percentile mode')
        parser.add_argument('--percentile', default=0.8, type=float, help='expression threshold, 0.9 as 90% in percentile mode')
        # parser.add_argument('--output-all', action='store_true')
        parser.add_argument('--open', action='store_true')
        args, cmds = parser.parse_known_args(arguments)
    else:
        args = arguments

    logger = tkutil.instanciate_standard_logger(sys._getframe().f_code.co_name)
    verbose = args.verbose
    auto_open = args.open
    seurat_dir = args.seurat_dir
    if verbose:
        logger.setLevel(logging.DEBUG)

    chart_type = args.dim
    # fn_input = list(sorted(args.u))[0]#, key=lambda f:os.))[0]

    outdir = args.o
    mode2d = (chart_type == 2)
    percentile = args.percentile
    # percentile_mode = args.percentile
    fn_expression = args.e or os.path.join(seurat_dir, 'seurat.normalized.tsv')
    fn_coord = args.u or os.path.join(seurat_dir, 'seurat.umap.tsv')
    fn_cluster = args.c or os.path.join(seurat_dir, 'seurat.clusters.tsv')
    # output_all = args.output_all
    threshold = args.percentile
    genes = args.g

    # percentile_mode = threshold < 1


    os.makedirs(outdir, exist_ok=True)

    fn_info = os.path.join(outdir, 'run.marker.info')
    fn_scatter = os.path.join(outdir, 'markerexpression.html')
    fn_violin = os.path.join(outdir, 'violin.html')
    fn_swarm = os.path.join(outdir, 'swarmplot.html')
    fn_total = os.path.join(outdir, 'clusters.html')

    logger.info('loading coordinates from {}'.format(fn_coord))
    coord = pd.read_csv(fn_coord, sep='\t', index_col=0)

    # logger.info('loading coordinates')
    # coord = pd.read_csv(fn_coord, sep='\t', index_col=0)

    clusterdf = pd.read_csv(fn_cluster, sep='\t', index_col=0)
    n_clusters = np.max(clusterdf) + 1
    n_cells = clusterdf.shape[0]
    logger.info(f'{n_clusters} clusters loaded')
            
    markers = []
    for gene in genes:
        if os.path.exists(gene):
            with open(gene) as fi:
                for line in fi:
                    items = re.split('[\\s,]+')
                    if len(items[0]) > 0:
                        markers.append(items[0])
        else:
            markers.append(gene)
    logger.info('loading {} gene expression from {}\n'.format(len(markers), fn_expression))
    expr = load_gene_expression(fn_expression, markers, verbose=verbose, logger=logger)
    if expr.shape[1] != n_cells:
        raise Exception('loaded expression was not compatible clusters having different cells')
    marker_genes = sorted([g for g in markers if g in expr.index])
    logger.info('{} viable genes loaded'.format(len(marker_genes)))

    if verbose:
        logger.info('Expression percentile [0,50,90,95,99,100]')
        for gene in sorted(marker_genes):
            if gene not in expr.index: continue
            values = expr.loc[gene].values
            pt = np.percentile(values, [0, 50, 90, 95, 99, 100])
            logger.info('{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}'.format(gene, pt[0], pt[1], pt[2], pt[3], pt[4], pt[5]))
    marker_size = max(2, 8 - np.log10(n_cells))
    traces = []
    # put all cells
    if mode2d:
        traces.append(go.Scattergl(x=coord.values[:,0], y=coord.values[:,1], name=f'All_{n_cells}', mode='markers',
            marker=dict(size=2, color='lightgray')))
    else:
        traces.append(go.Scatter3d(x=coord.values[:,0], y=coord.values[:,1], z=coord.values[:,2], name=f'All_{n_cells}', mode='markers',
            marker=dict(size=2, color='lightgray')))
    for gene in marker_genes:
        values = expr.loc[gene].values.reshape(-1)
        thr = np.percentile(values, percentile * 100)
        idx = [i for i in range(n_cells) if values[i] > thr]
        # print(gene, percentile, thr)
        if len(idx) > 0:
            xy = coord.iloc[idx,:].values
            logger.info('{}\t{} cells\t{}'.format(gene, thr, str(np.mean(xy, axis=0))))
            if mode2d:
                traces.append(go.Scattergl(x=xy[:,0], y=xy[:,1], name=gene, mode='markers',
                    marker=dict(size=marker_size)))
            else:
                traces.append(go.Scatter3d(x=xy[:,0], y=xy[:,1], z=xy[:,2], name=gene, mode='markers',
                    marker=dict(size=marker_size)))
    fig = go.Figure(traces)
    fig.update_layout(plot_bgcolor='white')
    fig.update_xaxes(linewidth=2, linecolor='black', showgrid=False, showline=True, mirror='ticks', title=coord.columns[0])
    fig.update_yaxes(linewidth=2, linecolor='black', showgrid=False, showline=True, mirror='ticks', title=coord.columns[1])

    plotly.offline.plot(fig, filename=fn_scatter, auto_open=auto_open)
    if 1:
        return
    

    cluster2index = {}
    with open(fn_info, 'w') as fo:
        fo.write('input:{}\n'.format(fn_expression))
        fo.write('n_cells:{}\n'.format(coord.shape[0]))
        if clusters is not None:
            fo.write('clusers:{}\n'.format(n_clusters))
            for cn in set(clusters):
                index = [j for j in range(len(clusters)) if clusters[j] == cn]
                n_cells_cluster = len(index)
                if n_cells_cluster > 0:
                    cluster2index[cn] = index
                    fo.write('cluster_C{}={} cells\n'.format(cn, n_cells_cluster))
        if verbose:
            for i, marker in enumerate(marker_genes):
                fo.write('marker{}:{}\t{}\n'.format(i + 1, marker, marker in expr.index))

    n_charts = min(expr.shape[0], 15)
    n_rows = min(4, n_charts)# // 4
    n_cols = (n_charts + n_rows - 1) // n_rows
    symbols = ['circle', 'diamond', 'square', 'triangle-up', 'circle-open-dot', 
    'diamond-open-dot', 'square-open-dot', 'cross', 'triangle-left', 'triangle-left-open']
    symbols = ['circle', ]

    cvals = clusterdf.values[:,0].values.reshape(-1)
    for c in range(-1, n_clusters):
        if c < 0:
            idx = [i for i in range(n_cells) if cvals[i] < 0]
        else:
            idx = [i for i in range(n_cells) if cvals[i] == c]
        xy = coord[idx,:]


    if clusters is not None: #
        import plotly.express as px
        # Violin chart
        n_charts = len(marker_genes)
        # fig = plotly.subplots.make_subplots(rows=n_rows, cols=n_cols, subplot_titles=marker_genes)
        # fig.print_grid()
        traces = []
        n_clusters = max(clusters) + 1
        pxtable = []
        for i, g in enumerate(marker_genes):
            values = expr.loc[g].values
            row = (i // n_cols) + 1
            col = (i % n_cols) + 1
            if 0 < row <= n_rows and 0 < col <= n_cols:
                labels = ['C{}'.format(j+1) for j in range(n_clusters)]
                for j, c in enumerate(clusters):
                    if j >= len(values):
                        sys.stderr.write('{} is out of bounds from {}\n'.format(j, len(values)))
                        continue
                    if values[j] > 0:
                        pxtable.append([g, 'C{}'.format(c+1), values[j]])

        pxdf = pd.DataFrame(sorted(pxtable, key=lambda r_:'{:03d}:{}'.format(int(r_[1][1:]), r_[0])), columns=['Gene', 'Cluster', 'Count'])
        f = px.strip(data_frame=pxdf, x='Cluster', y='Count', log_y=True, labels='Gene', facet_col='Gene', facet_col_wrap=n_rows)
        f.update_layout(plot_bgcolor='white')
        f.update_xaxes(linewidth=2, linecolor='black', showgrid=False, showline=True, mirror='ticks')
        f.update_yaxes(linewidth=2, linecolor='black', showgrid=False, showline=True, mirror='ticks')
        plotly.offline.plot(f, filename=fn_swarm, auto_open=auto_open)
    
        # fig.update_xaxes(title='Cluster', linewidth=1, linecolor='black', showgrid=False, showline=True, mirror='ticks')
        # fig.update_yaxes(title='Expression', linewidth=1, linecolor='black', showgrid=False, showline=True, mirror='ticks')
        # fig.update_layout(plot_bgcolor='white')
        # plotly.offline.plot(fig, filename=fn_violin, auto_open=False)
        # cluseter chart
        # color chart 
        # ES red
        # TS green
        # PrES blue
        comin = np.min(coord, axis=0)
        comax = np.max(coord, axis=0)
        xmin, xmax, ymin, ymax = comin[0], comax[0], comin[1], comax[1]

        if output_all:
            traces = []
            i_ = 0
            if mode2d:
                layout = go.Layout(
                    xaxis=dict(title=coord.columns[0], linewidth=2, linecolor='black', showgrid=False, showline=True, mirror='ticks'), 
                    yaxis=dict(title=coord.columns[1], linewidth=2, linecolor='black', showgrid=False, showline=True, mirror='ticks'),
                    plot_bgcolor='white'            
                )
            else:
                layout = go.Layout()
            for cn in sorted(cluster2index.keys()):
                index = cluster2index[cn]
                name = 'C{} ({})'.format(cn, len(index))
                xyz = coord.values[index,:]
                if cn < 0:
                    name = 'Unclassified'
                    color = 'lightgray'
                else:
                    color = None
                texts = expr.columns[index]
                if mode2d:
                    trace = go.Scattergl(x=xyz[:,0], y=xyz[:,1], name=name, mode='markers',
                        text=texts,
                        marker=dict(size=4, color=color,
                        symbol=symbols[i_ % len(symbols)]
                        ))
                else:
                    trace = go.Scatter3d(x=xyz[:,0], y=xyz[:,1], z=xyz[:,2], 
                        name=name, mode='markers', text=texts,
                        marker=dict(size=4, color=color, symbol=symbols[i_ % len(symbols)])
                    )
                traces.append(trace)
                i_ += 1
            fig = go.Figure(traces, layout)
            plotly.offline.plot(fig, fn_total)

    traces = []
    xyz = coord.values

    if mode2d:
        # print(xyz.shape)
        traces.append(go.Scattergl(x=xyz[:,0], y=xyz[:,1], mode='markers',text=list(coord.index),
            marker=dict(size=3, color='lightgray'), name='{} cells'.format(coord.shape[0])))
            
        layout = go.Layout(
            xaxis=dict(title=coord.columns[0], linewidth=2, linecolor='black', showgrid=False, showline=True, mirror='ticks'), 
            yaxis=dict(title=coord.columns[1], linewidth=2, linecolor='black', showgrid=False, showline=True, mirror='ticks'),
            plot_bgcolor='white'            
        )
    else:
        layout = go.Layout(
            scene=dict(
                xaxis=dict(title=coord.columns[0]), 
                yaxis=dict(title=coord.columns[1]),
                zaxis=dict(title=coord.columns[2])
            )
        )
        raise Exception('no 3d')


    for i, marker in enumerate(marker_genes):
        if marker not in expr.index: continue
        values = expr.loc[marker].values
        if percentile_mode:
            pt = np.percentile(values, [0, threshold * 100])
            thr_value = pt[1]
        else:
            thr_value = threshold
        index = []
        for j, v in enumerate(values):
            if v > thr_value:
                index.append(j)
        if len(index) == 0:
            continue
        if verbose and percentile_mode:
            # print(pt, np.percentile(values, [0,10,25,50,75,90,100]))
            sys.stderr.write('{}\t>{:.1f}\t{}/{}\n'.format(marker, thr_value, len(index), len(values)))
        xyz = coord.values[index]
        texts = coord.index[index]
        if mode2d:
            # print(texts[0:10])
            traces.append(go.Scattergl(x=xyz[:,0], y=xyz[:,1], mode='markers',
            name=marker,text=list(texts),
            marker=dict(size=4, symbol=symbols[i % len(symbols)])))
    fig = go.Figure(traces, layout)
    plotly.offline.plot(fig, filename=fn_scatter, auto_open=auto_open)

    return {'scatter.html':fn_scatter,
        'info':fn_info,
        'violin.html':fn_violin,
        'total.html':fn_total   
    }

def display_cluster_map(arguments=None):
    """Show calculated scatter plot in HTML.
Arguments:
    --seurat-dir seurat.R output directory
    -c cluster file in tsv (default seurat_dir/seurat.clusters.tsv)
    -u coordinate file in tsv (default seurat_dir/seurat.umap.tsv)
    --by-prefix [separator]
    """
    # if arguments is None:
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', help='expression tsv')
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('-g', nargs='+')
    parser.add_argument('-c', help='cluster')
    parser.add_argument('--color', metavar='tsv file', help='Color map of each cluster')
    parser.add_argument('-u', nargs='+', help='coordinates such as UMAP or tSNE')
    parser.add_argument('--dim', type=int, default=2, choices=[2, 3], help='2D or 3D chart')
    parser.add_argument('-o', default='exprchart')
    parser.add_argument('-t', default=1, type=float, help='expression threshold')
    # parser.add_argument('--percentile', action='store_true')
    parser.add_argument('--seurat-dir', default='seurat')
    # parser.add_argument('--colormode', type=int, default=-1, help='color mode')
    parser.add_argument('--preset-color', action='store_true')
    parser.add_argument('--min-cluster-size', default=10, type=int)
    parser.add_argument('--prefix', default=None)
    parser.add_argument('--open', action='store_true')
    args, cmds = parser.parse_known_args()#arguments)
    # else:
    #     args = arguments
    
    verbose = args.verbose
    chart_type = args.dim
    outdir = args.o
    seurat_dir = args.seurat_dir

    auto_open = args.open
    fn_out = os.path.join(outdir, 'chart.clusters.html')
    if args.prefix is not None:
        fn_group = os.path.join(outdir, 'chart.prefix.html')
        separator = args.prefix
    else:
        separator = None

    os.makedirs(outdir, exist_ok=True)
    # lower_limit = args.min_cluster_size

    fn_cluster = args.c or os.path.join(seurat_dir, 'seurat.clusters.tsv')
    fn_coord = args.c or os.path.join(seurat_dir, 'seurat.umap.tsv')
    fn_expression = args.e or os.path.join(seurat_dir, 'seurat.count.tsv')
    logger = tkutil.get_logger(sys._getframe().f_code.co_name)
    if args.verbose: logger.setLevel(10)
    # fn_log = os.path.join(outdir, 'chart.cluster')

    fn_info = os.path.join(outdir, 'run.info')
    fn_scatter = os.path.join(outdir, f'scatter.{chart_type}.html')
    fn_violin = os.path.join(outdir, 'violin.html')
    fn_total = os.path.join(outdir, 'clusters.html')

    info = {
        'scatter':fn_scatter,
        # 'violin':fn_violin,
        'total':fn_total,
    }
    print(chart_type)

    clu = pd.read_csv(fn_cluster, sep='\t', index_col=0)
    coord = pd.read_csv(fn_coord, sep='\t', index_col=0)
    if clu.shape[0] != coord.shape[0]:
        raise Exception('incompatible data size, number of cells of cluster labels and that of coordinates should match')
    n_cells = clu.shape[0]
    cluster2color = None
    n_clusters = np.max(clu.values) + 1
    min_clus = np.min(clu.values)

    logger.info('Cluster {}-{}'.format(min_clus, n_clusters))

    traces = []
    clu_vals = clu.iloc[:,0].values.reshape(-1)
    # marker_size = max(1, 10 - (np.log10(n_cells)))
    marker_size = min(6, max(1, 200 / (1 + np.sqrt(n_cells))))
    logger.info('{} cells, size={}'.format(n_cells, marker_size))
    for cn in range(min_clus, n_clusters):
        idx = [i for i in range(n_cells) if cn == clu_vals[i]]
        if cn < 0:
            cluster_name = 'U{}'.format(-cn)
        else:
            cluster_name = 'C{}_{}'.format(cn + 1, len(idx))
        logger.info(cluster_name)
        xy = coord.values[idx,:]
        x = xy[:,0]
        y = xy[:,1]
        if chart_type == 2:
            traces.append(go.Scattergl(x=x, y=y, name=cluster_name, mode='markers', marker=dict(size=marker_size)))
        else:
            z = xy[:,2]
            traces.append(go.Scatter3d(x=x, y=y, z=z, name=cluster_name, mode='markers', marker=dict(size=marker_size)))
    layout = go.Layout(plot_bgcolor='white')
    fig = go.Figure(traces, layout)        
    fig.update_xaxes(linewidth=2, linecolor='black', showgrid=False, showline=True, mirror='ticks', title=coord.columns[0])    
    fig.update_yaxes(linewidth=2, linecolor='black', showgrid=False, showline=True, mirror='ticks', title=coord.columns[1])    

    plotly.offline.plot(fig, filename=fn_scatter, auto_open=auto_open)

    if separator is not None:
        groups = {}
        for i, b in enumerate(clu.index):
            p_ = b.split(separator)[0]
            if p_ not in groups:
                groups[p_] = []
            groups[p_].append(i)

        traces = []
        clu_vals = clu.iloc[:,0].values.reshape(-1)
        marker_size = min(6, max(1, 100 / (1 + np.sqrt(n_cells))))
        logger.info('{} cells, size={}'.format(n_cells, marker_size))
        for group, idx in groups.items():
            # idx = [i for i in range(n_cells) if cn == clu_vals[i]]
            cluster_name = '{}_{}'.format(group, len(idx))
            logger.info(cluster_name)
            xy = coord.values[idx,:]
            x = xy[:,0]
            y = xy[:,1]
            if mode2d:
                traces.append(go.Scattergl(x=x, y=y, name=cluster_name, mode='markers', marker=dict(size=marker_size)))
            else:
                z = xy[:,2]
                traces.append(go.Scatter3D(x=x, y=y, z=z, name=cluster_name, mode='markers', marker=dict(size=marker_size)))
        layout = go.Layout(plot_bgcolor='white')
        fig = go.Figure(traces, layout)        
        fig.update_xaxes(linewidth=2, linecolor='black', showgrid=False, showline=True, mirror='ticks', title=coord.columns[0])    
        fig.update_yaxes(linewidth=2, linecolor='black', showgrid=False, showline=True, mirror='ticks', title=coord.columns[1])    

        plotly.offline.plot(fig, filename=fn_group, auto_open=auto_open)


def display_count_map(arguments=None):
    """Draw count distribution charts by cluster

    return
    outputs = {
        'count.tsv':fn_count,
        'stats.tsv':fn_stat,
        'color.tsv':fn_color,
        'depth.graph.html':fn_depth,
        'count.graph.html':fn_histogram
    }
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', metavar='TSV file', help='cluster TSV file')
    parser.add_argument('-m', '--mtx', metavar='directory', help='matrix directory')
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--forced', action='store_true')
    parser.add_argument('-o', default='stat')
    parser.add_argument('-e', default=None)
    parser.add_argument('-u', default=None)    
    # parser.add_argument('--bin', type=int, default=50)
    parser.add_argument('--seurat-dir')
    parser.add_argument('--open', action='store_true')
    parser.add_argument('--barcode-filter', default=None, help='prefix of using barcodes')
    parser.add_argument('--mitochondria', type=float, default=5.0, help='threshold mitochondria RNA (%)')
    if arguments is not None:
        args = parser.parse_known_args(arguments)[0]
    else:
        args = parser.parse_knwon_args()[0]
    basedir = args.seurat_dir
    outdir = args.o or basedir
    fn_count = args.e or os.path.join(basedir, 'seurat.count.tsv')
    fn_cluster = args.c or os.path.join(basedir, 'seurat.clusters.tsv')
    fn_coord = args.u or os.path.join(basedir, 'seurat.umap.tsv')
    mt_thr = args.mitochondria
    mtdir = args.mtx
    logger = tkutil.instanciate_standard_logger(sys._getframe().f_code.co_name)
    auto_open = args.open

    if args.verbose: 
        logger.setLevel(10)

    def serialize_barcodes(bcd):
        return re.sub('[^A-Za-z]', '', bcd.strip('"'))

    logger.info('loadin clusters')
    df_clus = pd.read_csv(fn_cluster, sep='\t', index_col=0)
    df_clus.index = [serialize_barcodes(b) for b in df_clus.index]

    fn_cache = os.path.join(basedir, '.count.cache')
    logger.info(f'count cache is {fn_cache}')
    columns = ['n_reads', 'n_features', 'n_mt']
    count_df = None
    if os.path.exists(fn_cache):
        count_df = pd.read_csv(fn_cache, sep='\t', index_col=0)
        if len(count_df.columns) != len(columns) or count_df.shape[0] != df_clus.shape[0]:
            count_df = None
    if count_df is None: # read
        if mtdir is not None and os.path.exists(mtdir):
            logger.info('loading counts from sparse matrix')
            mtx = scipy.io.mmread(mtdir).tocsr()
            n_reads = np.asarray(mtx.sum(axis=0)).reshape(-1)
            features = cellrangerwrapper.load_features(mtdir)
            barcodes = [serialize_barcodes(b) for b in cellrangerwrapper.load_barcodes(mtdir)]
            midx = []
            for i, f in enumerate(features):
                if f.lower().startswith('mt-'):
                    midx.append(i)
            n_mt = np.asarray(mtx[midx,:].sum(axis=0)).reshape(-1)
            mtx[mtx>0] = 1
            n_ft = np.asarray(mtx.sum(axis=0)).reshape(-1)
        else:
            logger.info('loading counts from csv')
            df = pd.read_csv(fn_count, sep='\t', index_col=0)
            features = df.index
            barcodes = [serialize_barcodes(b) for b in df.columns]
            midx = []
            for i, f in enumerate(features):
                if f.lower().startswith('mt-'):
                    midx.append(i)
            n_reads = np.sum(df.values, axis=0)
            n_mt = np.sum(df.iloc[midx], axis=0)
            df[df>0] = 1
            n_features = np.sum(df, axis=0)
            count_df = pd.DataFrame([n_reads, n_features, n_mt], columns=barcodes, index=columns).T
        logger.info(f'saving to {fn_cache}')
        count_df.to_csv(fn_cache, sep='\t')
    clusters = set(df_clus.values.reshape(-1))
    n_clusters = max(clusters) + 1
    count_dist = {}
    mt_dist = {}
    idx = {}
    for c in clusters:
        if c not in count_dist: 
            count_dist[c] = []
            mt_dist[c] = []
            idx[c] = []
    for i, c in enumerate(df_clus.values[:,0].reshape(-1)):
        idx[c].append(i)

    # MT and read count map
    fn_chart = os.path.join(outdir, 'count_mt.html')
    figs = plotly.subplots.make_subplots(rows=2, cols=2, subplot_titles=['Count', 'MT', 'Clusters'])
    mtratio = np.asarray(count_df.iloc[:,2]  / count_df.iloc[:,0]).reshape(-1)
    tot_avg = np.mean(count_df.iloc[:,0])
    tot_sd = np.std(count_df.iloc[:,0], axis=0)
    mt_avg = np.mean(mtratio)
    mt_sd = np.std(mtratio)
    logger.info('Total avg={:.0f}, sd={:.0f}, MT avg={:.1f}%, sd={:.2f}%'.format(
        tot_avg, tot_sd, mt_avg * 100, mt_sd * 100
    ))
    coord = pd.read_csv(fn_coord, sep='\t', index_col=0)
    fn_stat = os.path.join(outdir, 'readscounts.tsv')
    fo = open(fn_stat, 'w')
    fo.write('\tn_cells\tn_reads\tn_features\tn_mt\tpercent_mt\n')
    
    for c in sorted(clusters):
        nr = np.asarray(count_df.iloc[idx[c],0].values).reshape(-1)
        ft_mean = np.mean(count_df.iloc[idx[c], 1])

        cnt_mean = np.mean(nr)
        cnt_z = (cnt_mean - tot_avg) / tot_sd

        mt_ratios = mtratio[idx[c]]
        mt_mean = np.mean(mt_ratios)
        mt_z = (mt_mean - mt_avg) / mt_sd

        cname = 'C{}'.format(c+1) if c >= 0 else 'unc'

        fo.write('{}\t{}\t{}\t{}\n'.format(cname, len(idx[c]), cnt_mean, ft_mean, mt_mean))

        ostr = '{}\t{}\t{:.0f}\t{:.2f}\t{:.0f}\t{:.2f}'.format(cname, len(idx[c]), cnt_mean, cnt_z, mt_mean * 100, mt_z)
        logger.info(ostr)
        xy = coord.values[idx[c],:]

        def __zcolor(z, zmin=-2, zmax=2, mode=0):
            if z < 0:
                deg = min(255, int((z / zmin) * 256))
                b = 224 + deg // 4
                r = 255 - deg
                g = max(0, 255 - deg * 2)
            else:
                deg = min(255, int((z / zmax) * 256))
                r = 224 + deg // 4
                b = 255 - deg
                g = max(0, 255 - deg * 2)
            if mode == 0:
                return 'rgb({},{},{})'.format(r, g, b)
            else:
                return 'rgb({},{},{})'.format(b, r, g)
        figs.add_trace(go.Scattergl(x=xy[:,0], y=xy[:,1], name='{} {:.0f}'.format(cname, cnt_mean), mode='markers', 
            marker=dict(size=4, color=__zcolor(cnt_z)), legendgroup='1'), row=1, col=1)
        figs.add_trace(go.Scattergl(x=xy[:,0], y=xy[:,1], name='{} {:.1f}%'.format(cname, mt_mean * 100), mode='markers', 
            marker=dict(size=4, color=__zcolor(mt_z, mode=1)), legendgroup='2'), row=1, col=2)
        figs.add_trace(go.Scattergl(x=xy[:,0], y=xy[:,1], 
            name='{} {} cells'.format(cname, len(idx[c])), mode='markers', 
            marker=dict(size=4), legendgroup='3'), row=2, col=1)
    figs.update_xaxes(linewidth=2, linecolor='black', showgrid=False, showline=True, mirror='ticks', title=coord.columns[0])
    figs.update_yaxes(linewidth=2, linecolor='black', showgrid=False, showline=True, mirror='ticks', title=coord.columns[1])
    figs.update_layout(plot_bgcolor='white', legend_tracegroupgap = 100)
    # figs.update_layout(plot_bgcolor='white', showlegend=True, row=2, col=1)
    plotly.offline.plot(figs, filename=fn_chart, auto_open=auto_open)
        

def count_tags_by_cluster(arguments=None):
    """Draw count distribution charts by cluster

    return
    outputs = {
        'count.tsv':fn_count,
        'stats.tsv':fn_stat,
        'color.tsv':fn_color,
        'depth.graph.html':fn_depth,
        'count.graph.html':fn_histogram
    }
    """
    raise Exception('buggy')
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', metavar='TSV file', help='cluster TSV file')
    parser.add_argument('-m', '--mtx', metavar='directory', help='matrix directory')
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--forced', action='store_true')
    parser.add_argument('-o', default='stat')
    parser.add_argument('-e', default=None)
    parser.add_argument('--bin', type=int, default=50)
    parser.add_argument('--seurat-dir')
    parser.add_argument('--open', action='store_true')
    parser.add_argument('--barcode-filter', default=None, help='prefix of using barcodes')
    parser.add_argument('--mitochondria', type=float, default=5.0, help='threshold mitochondria RNA (%)')
    if arguments is not None:
        args = parser.parse_known_args(arguments)[0]
    else:
        args = parser.parse_knwon_args()[0]
    # print(arguments)

    verbose = args.verbose
    seurat_dir = args.seurat_dir
    auto_open = args.open
    matrix_dir = args.mtx or os.path.join(os.path.dirname(seurat_dir), 'mtx')
    logger = tkutil.instanciate_standard_logger(sys._getframe().f_code.co_name)
    prefix = args.barcode_filter
    mit_thr = args.mitochondria
    if args.verbose: logger.setLevel(logging.DEBUG)

    # if srcdir is None:
    #     raise Exception('no sparse directory given with --mtx option')

    # set filenames of output
    outdir = args.o
    os.makedirs(outdir, exist_ok=True)
    fn_coord = args.c or os.path.join(outdir, 'seurat.umap.tsv')
    fn_count = os.path.join(outdir, 'count_by_cluster.tsv')
    fn_stat = os.path.join(outdir, 'count_stat.tsv')
    fn_color = os.path.join(outdir, 'cluster_color.tsv')
    fn_depth = os.path.join(outdir, 'depth.html')
    fn_histogram = os.path.join(outdir, 'count_by_cluster.html')
    fn_expr = args.e or os.path.join(seurat_dir, 'seurat.count.tsv')
    fn_cluster = args.c or os.path.join(seurat_dir, 'seurat.clusters.tsv')
    fn_map_chart = os.path.join(outdir, 'count_by_cluster.html')
    forced = args.forced

    if os.path.exists(fn_stat) is False or forced:
        if verbose: sys.stderr.write('\033[Kloading clusters\r')
        clusterdf = pd.read_csv(fn_cluster, sep='\t', index_col=0)
        clusterdf.index = [re.sub('[^A-Z]', '', bcd) for bcd in clusterdf.index]
        # print(clusterdf)
        if 'Cluster' in clusterdf.columns:
            clusters = clusterdf['Cluster'].values
        else:
            clusters = clusterdf.values[:,0]
        n_clusters = max(clusters) + 1
        n_cells = clusterdf.shape[0]
        fn_out = args.o
        logger.info('{} clusters of {} cells loaded'.format(n_clusters, n_cells))

        if os.path.exists(matrix_dir) is False:
            logger.info(f'no matrix directory count from table from {fn_expr}')
            # cdf = load_reads_from_count_matrix(fn_expr)
            cdf = pd.read_csv(fn_expr, sep='\t', index_col=0)

            # def __reformat_barcode(barcodes):
            #     bcd = []
            #     for b in barcodes:
            #         bcd.append(re.sub('(DK\\d+_\\d+)\\.([ACGT]+)\\.(\\d+)', '\\1:\\2-\\3', b.strip('"')))
            #         pass
            #     return bcd
            cdf.index = [re.sub('[^A-Z]', '', bcd) for bcd in cdf.index]#__reformat_barcode(cdf.index)
            barcodes = cdf.index
            # cdf = load_total_counts(fn_expr)
        else:
            logger.info(f'loading barcodes from {matrix_dir}')
            barcodes = load_barcodes(matrix_dir)
        barcode2column = {}
        for i, b in enumerate(barcodes):
            b = re.sub('[^A-Z]', '', b.strip('"'))
            barcode2column[b] = i
        # print(list(barcode2column.keys())[0:10])
        cluster2bcd = {}
        for cn in range(n_clusters):
            cluster2bcd[cn] = []    
        for i, barcode in enumerate(clusterdf.index):
            if clusters[i] >= 0:
                try:
                    cluster2bcd[clusters[i]].append(barcode2column[barcode])
                except:
                    print(barcode, barcode in barcode2column)
                    raise
        if verbose:
            sys.stderr.write('\033[K{} barcodes loaded\n'.format(len(barcodes)))

        if verbose:
            sys.stderr.write('\033[Kloading matrix\r')
        if cdf is None:
            cdf = load_reads_from_sparse_matrix(matrix_dir)
        counts_per_cell = cdf['n_Reads'].values.reshape(-1)
        features_per_cell = cdf['n_Features'].values.reshape(-1)
        mt_per_cell = cdf['n_MT'].values.reshape(-1)
        if verbose:
            sys.stderr.write('\033[K{} x {} matrix loaded\n'.format(cdf.shape[0], cdf.shape[1]))

        n_cells = len(barcodes)
        reads_and_genes = np.zeros((n_cells, 3), dtype=np.float32)
        for i, barcode in enumerate(barcodes):
            col = barcode2column[barcode]
            if barcode in clusterdf.index:
                clus = clusterdf.loc[barcode].values[0]
            else:
                clus = -1
            n_reads = counts_per_cell[i]
            n_features = features_per_cell[i]
            n_mt = mt_per_cell[i]
            reads_and_genes[i] = (clus, n_reads, n_features)
        stat = pd.DataFrame(reads_and_genes, columns=['Cluster', 'n_reads', 'n_features', 'n_MT'], index=barcodes, dtype=np.int32)
        stat.to_csv(fn_stat, sep='\t')
    else:
        stat = pd.read_csv(fn_stat, sep='\t', index_col=0)

    # filter barcodes
    if prefix is not None:
        idx = []
        for i, b in enumerate(stat.index):
            if b.startswith(prefix):
                idx.append(i)
        logger.info('{} in {} barcodes were selected'.format(len(idx), stat.shape[0]))
        stat = stat.iloc[idx, :]

    # scatter by cluster
    import collections
    clusters = set(stat['Cluster'])
    cluster2index = {}
    for i, c in enumerate(stat['Cluster'].values):
        if c not in cluster2index: cluster2index[c] = []
        cluster2index[c].append(i)
    cluster_to_data = collections.OrderedDict()

    fn_scatter = os.path.join(outdir, 'coloredclusters.html')
    traces = []
    cluster_title = []

    field_of_interest = 'n_reads' # color
    total = stat[field_of_interest].values
    total_mean = np.mean(total)
    total_sd = np.std(total)
    total_med = np.median(total)
    cluster_color = {}
    print(stat)
    if mit_thr > 0:
        mt = stat['n_MT'].values.reshape(-1)
        tot = stat['n_reads'].values.reshape(-1)
        mt_ratio = mt * 100 / tot
    else:
        mt_ratio = None

    for cn in sorted(clusters):
        submat = stat[stat['Cluster']==cn]
        cluster_to_data[cn] = submat
        X = submat['n_reads']
        Y = submat['n_features']

        values = submat[field_of_interest].values
        submat_mean = np.mean(values)
        submat_sd = np.std(values)
        submat_med = np.median(values)
        z = (submat_mean - total_mean) / total_sd
        deg = z
        u = min(255, max(0, int((deg + 1) * 256)))
        d = min(255, max(0, int((1 - deg) * 256)))
        red = 255 if deg > 0 else u
        green = min(u, d)
        blue = 255 if deg < 0 else d
        if mt_ratio is not None:
            mt_med = np.median(mt_ratio[cluster2index[cn]])
            green = min(1, (mt_med / mit_thr)) * 64 + 128

        cluster_color[cn] = ['rgb({:.0f},{:.0f},{:.0f})'.format(red, green, blue), z]
        if cn >= 0:
            title = 'C{} ({})'.format(cn + 1, submat.shape[0])
            cluster_title.append(title)
        else:
            title = 'Unclassified ({})'.format(submat.shape[0])
        traces.append(go.Scatter(x=X, y=Y, text=submat.index,
            mode='markers', 
            marker=dict(size=5),
            name=title))
    # set color
    if cluster_color is not None and len(cluster_color) > 0:
        with open(fn_color, 'w') as fo:
            for cn in sorted(clusters):
                cinfo = cluster_color[cn]
                fo.write('{:.0f}\t{}\t{:.4f}\n'.format(cn, cinfo[0], cinfo[1]))

    # counts vs features
    fig = go.Figure(traces)
    fig.update_layout(plot_bgcolor='white')
    fig.update_xaxes(title='n_counts', linewidth=2, linecolor='black', showgrid=False, showline=True, mirror='ticks')
    fig.update_yaxes(title='n_genes', linewidth=2, linecolor='black', showgrid=False, showline=True, mirror='ticks')
    plotly.offline.plot(fig, filename=fn_depth, auto_open=False)

    # cluster map
    #display_cluster(map)
    # coord = pd.read_csv(fn_coord, sep='\t', index_col=0)
    # traces = []
    # for cn in sorted(clusters):
    #     barcodes = stat[stat['Cluster']===cn].index
    #     xyz = coord.loc[barcodes]
    #     traces.append(go.Scattergl(x=xyz[:,0], y=xyz[:,1], name='C{} {}'.format(), mode='markers',
    #         marker=dict(size=5, color=cluster_color.get(cn, 'gray'))))

    # mitochondria
    # if mit_thr > 0:
    #     mt = stat['n_mt'].values.reshape(-1)
    #     tot = stat['n_reads'].values.reshape(-1)
    #     percent = mt / tot
    #     mcount = stat[['n_mt', 'n_reads']].values
    #     idx = []
    #     for i, b in enumerate(stat.index):
    #         vals = mcount[i,:]
    #         if vals[0] * mit_thr / 100 > vals[1]:
    #             idx.append(i)
    #     if len(idx) > 0:
    #         trace = go.Scatter(x=)
    #         fig.add_trace(trace, row=row, col=col)
        
    #     traces.append(go.Scatter(x=X, y=Y, name='MT >={:.1}%'.format(mit_thr), mode='markers', 
    #     marker=dict(symbol='ciecle-open', size=7, color='darkgreen')))

    # histogram
    binsize = args.bin
    n_clusters = len(cluster_title)#df.shape[1]
    n_chart_cols = 3
    n_chart_rows = (n_clusters + 1 + n_chart_cols - 1) // n_chart_cols
    maxcnt = ((np.max(stat['n_reads']) + binsize - 1) // binsize) * binsize

    fig = plotly.subplots.make_subplots(
        cols=n_chart_cols, 
        rows=n_chart_rows, 
        subplot_titles=cluster_title)

    xlimit = min(30000, maxcnt) #// binsize
    gsize = xlimit // binsize
    x = np.arange(0, maxcnt, binsize)
    accum = np.zeros(x.shape[0], dtype=np.int32)

    index = 0
    for cn in sorted(cluster_to_data.keys()):
        if cn < 0: continue
        row = (index // n_chart_cols) + 1
        col = (index % n_chart_cols) + 1
        y = np.zeros(x.size, dtype=np.int32)
        submat = cluster_to_data[cn]
        n_reads = submat['n_reads']
        for n in submat['n_reads'].values:
            y[n // binsize] += 1
        accum += y
        trace = go.Scatter(x=x[0:gsize], y=y[0:gsize], name=cluster_title[index], mode='lines', line=dict(color='navy'))
        fig.add_trace(trace, row=row, col=col)
        index += 1


    trace = go.Scatter(x=x, y=accum, name='Total')
    fig.add_trace(trace, row=n_chart_rows, col=n_chart_cols)

    fig.update_xaxes(range=[0,xlimit])
    fig.update_xaxes(showgrid=False, showline=True, linewidth=2, linecolor='black')
    fig.update_yaxes(showgrid=False, showline=True, linewidth=2, linecolor='black')
    
    fig.update_layout(legend=None, plot_bgcolor='white')
    plotly.offline.plot(fig, filename=fn_histogram, auto_open=auto_open)

    outputs = {
        'cluter.count.tsv':fn_count,
        'stats.tsv':fn_stat,
        'color.tsv':fn_color,
        'depth.graph.html':fn_depth,
        'count.graph.html':fn_histogram
    }

    return outputs

def generate_cluster_heatmap(arguments=None):
    if arguments is None:
        parser = argparse.ArgumentParser()
        parser.add_argument('--verbose', action='store_true')
        parser.add_argument('-s', '--sep', default='_', help='delimiter between name and barcode')
        parser.add_argument('-u', help='coordinates such as UMAP or tSNE')
        parser.add_argument('-o', default='chart')
        parser.add_argument('--open', action='store_true')
        args = parser.parse_known_args(arguments)[0]

        parser = argparse.ArgumentParser()
        parser.add_argument('-c', metavar='TSV file', help='cluster TSV file')
        parser.add_argument('-i', metavar='directory', help='matrix directory')
        parser.add_argument('--verbose', action='store_true')
        parser.add_argument('--forced', action='store_true')
        parser.add_argument('-o', default='stat')
        parser.add_argument('-e', default=None)

        parser.add_argument('--seurat-dir')
        parser.add_argument('-s', default=None, help='sample separator')
        parser.add_argument('-g', nargs='+')
        
        parser.add_argument('--open', action='store_true')

        args = parser.parse_known_args(arguments)[0]

        verbose = args.verbose
        srcdir = args.i
    else:
        args = arguments

    auto_open = args.open
    outdir = args.o
    verbose = args.verbose
    seurat_dir = args.seurat_dir
    # if seurat_dir is not None:
    #     fn_clus = 
    #     fn_expr = os.path.join(seurat_dir, 'seurat.normalized.tsv')

    logger = tkutil.instanciate_standard_logger(sys._getframe().f_code.co_name)
    if args.verbose: logger.setLevel(10)

    fn_clus = args.c or os.path.join(seurat_dir, 'seurat.clusters.tsv')
    fn_expr = args.e or os.path.join(seurat_dir, 'seurat.normalized.tsv')

    fn_heatmap = os.path.join(outdir, 'heatmap.html')

    clusters = pd.read_csv(fn_clus, sep='\t', index_col=0)
    barcodes = clusters.index
    genes = args.g
    n_clusters = np.max(clusters.values) + 1
    logger.info(f'{n_clusters} clusters found')
    outdir = args.o
    os.makedirs(outdir, exist_ok=1)

    separator = args.s
    import collections
    if separator is not None:
        samples = collections.OrderedDict()
        for i, b in enumerate(barcodes):
            if b.find(separator) >= 0:
                sample = b[:b.find(separator)]
                if sample not in samples: 
                    logger.info(f'new sample \t{sample}')
                    samples[sample] = [i, ]
                else:
                    samples[sample].append(i)
            else:
                raise Exception(f'barcode {b} does not have separator')
    else:
        samples = {'Sample':barcodes}

    expr = load_gene_expression(fn_expr, genes, verbose=verbose, logger=logger)
    # print(expr)
    # exit(0)
    genes = [g for g in genes if g in expr.index]

    # marker expression
    n_rows = len(genes)
    n_cols = n_clusters
    c2i = [[] for i in range(n_clusters)]
    # print(c2i)
    columns = []
    cells = []
    for i in range(clusters.shape[0]):
        c = clusters.values[i,0]
        # print(i, c)
        # logger.info(i, clusters.iloc[i,0], clusters.values[i,0])
        if c >= 0:
            c2i[c].append(i)
    c2e = np.array((n_rows, n_cols), dtype=np.float32) # cluster expression
    composition = np.zeros((len(samples), n_clusters), dtype=np.int32)
    for c in range(n_clusters):
        columns += c2i[c]
        if samples is None:
            cells += ['C{}_{}'.format(c + 1, i + 1) for i in range(c2i[c])]
        else:
            scnt = {}
            for i in c2i[c]:
                barcode = barcodes[i]
                if separator is None or barcode.find(separator) is None:
                    sn = ''
                else:
                    sn = barcode[:barcode.find(separator)]
                scnt[sn] = scnt.get(sn, 0) + 1
                cells.append(f'{sn}_C{c+1}_{scnt[sn]}')
            for i, s in enumerate(samples.keys()):
                composition[i,c] = len([j for j in samples[s] if j in c2i[c]])
            # cells += [expr.columns[i] for i in c2i]
    # print(composition)
    fn_composition = os.path.join(outdir, 'cluster_composition.tsv')
    df = pd.DataFrame(composition, index=samples.keys(), columns=[f'C{i+1}' for i in range(n_clusters)])
    df.to_csv(fn_composition, sep='\t')

    sort_index = []
    ctable = []
    y = []
    # score_matrix = []
    rank_matrix = np.zeros((len(genes), n_clusters), dtype=np.int)
    for i, gene in enumerate(genes):
        values = expr.loc[gene].values
        row = np.zeros(n_clusters, dtype=np.float)
        valset = []
        for c in range(n_clusters):
            idx = c2i[c]
            if len(idx) > 0:
                val = values[idx]
                row[c] = np.mean(val)
                valset.append(val)
        pval = scipy.stats.f_oneway(*valset).pvalue # ANOVA
        sval = list(sorted(row))
        maxval = sval[-1]
        repindex = sval[-1] - sval[-2]
        ctable.append(row)
        if pval > 0.01: 
            repindex = 0
            continue
        sort_index.append(repindex)
        y.append(gene)
        # sr = np.zeros(n_clusters)
        for rank, c in enumerate(sorted(range(n_clusters), key=lambda c_:row[c_], reverse=True)):
            if row[c] * 4 > maxval:
                rank_matrix[i, c] = n_clusters - rank
            else:
                break
            # sr[c] = rank
            # sr.append((len(clusters) - rank + 1) * (2 ** i))
        # score_matrix.append(sr)
    df = pd.DataFrame(rank_matrix, index=genes, columns=[f'C{i+1}' for i in range(n_clusters)])
    # print(df)
    # scores = np.arange(1, len(genes)+ 1).reshape(-1, 1) * df
    # print(scores)
    touched = set()
    sorted_clusters = []
    for i, g in enumerate(genes):
        sc = sorted(range(n_clusters), key=lambda c_:rank_matrix[i,c_], reverse=True)
        # print(i, g, sc)
        for c in list(sc)[0:6]:
            if rank_matrix[i,c] == 0:
                break
            if c not in touched:
                logger.info('{} th cluster {} for {}'.format(len(sorted_clusters), c, g))
                sorted_clusters.append(c)
                touched.add(c)
            

    # print(sorted_clusters)
    colindex = []
    cells = []
    
    for c in sorted_clusters:
        columns += c2i[c]
        # print('{} C{}'.format(len(cells), c+1))
        # print(c2i[c])

        if samples is None:
            cells += ['C{}_{}'.format(c + 1, i + 1) for i in range(len(c2i[c]))]
        else:
            scnt = {}
            for i in c2i[c]:
                barcode = barcodes[i]
                if separator is None or barcode.find(separator) is None:
                    sn = ''
                else:
                    sn = barcode[:barcode.find(separator)]
                scnt[sn] = scnt.get(sn, 0) + 1
                cells.append(f'{sn}_C{c+1}_{scnt[sn]}')
        idx = c2i[c]
        colindex += idx
    # print(len(colindex))
    expr = expr.iloc[:,colindex]
    expr.columns = cells

    z = []
    x = cells
    y = genes
    for i, g in enumerate(genes):
        z.append(expr.loc[g].values)#[columns])
    traces = []
    traces.append(go.Heatmap(z=z, x=x, y=y, zmin=0, colorscale='reds'))
    fig = go.Figure(traces)
    fig['layout']['yaxis']['autorange'] = 'reversed'
    plotly.offline.plot(fig, filename=fn_heatmap, auto_open=auto_open)
    fn_res = os.path.join(outdir, 'raw_hatmap.tsv')

    # cell numbers


    df = pd.DataFrame(z, index=y, columns=x)
    df.to_csv(fn_res, sep='\t', float_format='%.4f')



def main():
    parser = argparse.ArgumentParser()
    # parser.add_argument('cmd', choices=['marker', 'count', 'cluster', 'byname', 'all'])
    parser.add_argument('--seurat-dir', default=None, metavar='directory', help='output directory of seurat.R script')
    parser.add_argument('-e', help='expression tsv', default=None)
    parser.add_argument('-g', nargs='+', default=None, help='marker genes for marker command')
    parser.add_argument('-u', default=None, help='coordinates such as UMAP or tSNE')
    parser.add_argument('-o', default='analysis', metavar='directory', help='output')

    parser.add_argument('-t', metavar='number', default=1, type=float, help='expression threshold')
    parser.add_argument('--percentile', type=float, default=1.0)

    parser.add_argument('-c', metavar='filename', help='cluster TSV file')
    parser.add_argument('-i', metavar='directory', help='matrix directory')

    parser.add_argument('--forced', action='store_true', help='force calculation')
    parser.add_argument('--bin', type=int, default=100, metavar='number', help='count command graph bin')

    parser.add_argument('-s', default='_', metavar='char', help='Sample prefix separator')
    parser.add_argument('--color', metavar='tsv file', help='Color map of each cluster')

    parser.add_argument('--dim', type=int, default=2, choices=[2, 3], help='2D or 3D chart')
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--open', action='store_true')

    parser.add_argument('-m', '--mtx', metavar='directory', help='matrix directory')
    parser.add_argument('--mitchondria', type=float, default=5.0, help='threshold mitochondria RNA in percent')

    parser.add_argument('--dime', type=int, default=2, choices=[2, 3], help='2D or 3D chart')

    args, commands = parser.parse_known_args()
    cargs = list(sys.argv[2:])
    # cmd = args.cmd
    if args.seurat_dir is not None:
        basedir = args.seurat_dir
        # if args.e is None: 
        #     cargs += ['-e', os.path.join(basedir, 'seurat.normalized.tsv')]
        if args.u is None:
            cargs += ['-u', os.path.join(basedir, 'seurat.umap.tsv')]
        if args.c is None:
            cargs += ['-c', os.path.join(basedir, 'seurat.clusters.tsv')]

    for cmd in commands:
        if cmd == 'all':
            if args.g is not None:
                show_marker(cargs)
            results = count_tags_by_cluster(cargs)
            if args.e is None:
                cargs += ['-e', results['count.tsv']]
            if args.g is not None:
                show_marker(args)
            cargs += ['--color', results['color.tsv']]
            display_cluster_map(args)
        elif cmd == 'marker':
            show_marker(args)
        elif cmd == 'count':
            info = display_count_map(cargs)
        elif cmd == 'cluster':
            display_cluster_map(args)
        elif cmd == 'byname':
            display_cell_distribution(cargs)
        elif cmd == 'heatmap':
            generate_cluster_heatmap()
        elif cmd == 'swarm' or cmd == 'swarmplot':
            draw_swarmplot()

        # else:
        #     raise Exception('not implemented {}'.format(cmd))

if __name__ == '__main__':
    main()
