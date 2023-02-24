import os, sys, re, argparse
import logging
import functools
import json

sys.path.append('/home/data/cellranger/scripts/')

import tkutil

class AnnotatedGene(object):
    def __init__(self, line):
        if line == None:
            self.chromosome = None
            self.start = -1
            self.end = -1
            self.gene_name = None
            self.gene_number = -1
            self.strand = '.'
            self.coding_start = -1
            self.coding_end = -1
            self.exon_lengths = []
            self.exon_start = []
            self.exon_end = []
            self.total_length = 0
            self.biotype = None
            return
            
        items = line.strip().split('\t')
        self.chromosome = re.sub('^chr', '',  items[1])
        self.start = int(items[2])
        self.end = int(items[3])
        m = re.match('(.*)\\.(\\d+)$', items[4])
        if m:
            self.gene_name = m.group(1)
            self.transcript_id = items[4]
        else:
            self.gene_name = self.transcript_id = items[4]
        self.strand = items[6]
        self.coding_start = int(items[7])
        self.coding_end = int(items[8])
        n_exons = int(items[10])
        self.exon_lengths = [int(x_) for x_ in items[11].strip(',').split(',')]
        self.exon_start = [int(x_) for x_ in items[12].strip(',').split(',')]
        exon_positions = []
        pos = self.start
        for i in range(n_exons):
            exon_positions.append(self.exon_start[i] + self.exon_lengths[i])
        self.exon_end = exon_positions
        self.total_length = sum(self.exon_lengths)
        self.biotype = None

    def __is_coding(self):
        if self.biotype is not None:
            return self.biotype == 'protein_coding'
        return self.coding_start >= 0
    is_coding = property(__is_coding)
    def find_exon(self, chromosome, position):
        n_exons = len(self.exon_start)
        if self.chromosome != chromosome or position < self.start or self.end < position:
            return None
        for i in range(n_exons):
            if self.exon_start[i] <= position <= self.exon_end[i]:
                exon_num = i if self.strand == '+' else n_exons - 1 - i
                offset = sum(self.exon_lengths[:i]) + position - self.exon_start[i]
                if self.strand == '+':
                    offset = self.total_length - offset
                return [self, i, self.is_coding, offset]
        return None
    
    def __str__(self):
        return '{}\t{}:{}:{}-{}'.format(self.gene_name, self.chromosome, self.strand, self.start, self.end, len(self.exon_lengths), self.total_length)

    def __repr__(self):
        return '{}\t{}:{}:{}-{}\t{}\t{}'.format(self.gene_name, self.chromosome, self.strand, self.start, self.end, len(self.exon_lengths), self.total_length)
    
    @classmethod
    def load_gtf(cls, filename, **kwargs):
        """
"""
        import gzip, hashlib, pickle
        md5 = hashlib.md5()
        md5.update(os.path.abspath(filename).encode('utf-8'))
        logger = kwargs.get('logger', logging.getLogger())
        
        cache = os.path.join(os.path.dirname(filename), '.' + md5.hexdigest()[0:5] + '.cache')
        #forced = kwargs.get('forced', False)
        biotypes = kwargs.get('biotype', None)
        if isinstance(biotypes, str):
            biotypes = [biotypes, ]
        forced = kwargs.get('forced', False)
        #forced = True
        accepted_features = ['exon', 'start_codon', 'stop_codon']
        if tkutil.check_file(cache, filename) and not forced:
            logger.info('loading data from cache')
            with gzip.open(cache, 'rb') as fi:
                genes = pickle.load(fi)
                pass
            pass
        else:
            genes = {}
            n_tr = n_skip = 0
            gene_names = set()
            
            with open(filename) as fi:
                for line in fi:
                    items = line.strip().split('\t')
                    if len(items) < 9 or (items[2] not in accepted_features) or items[3].isdigit() is False:
                        n_skip += 1
                        continue
                    
                    #print(items)
                    elems = items[8].split(';')
                    prop = {}
                    for elem in elems:
                        m = re.match('\\s*(\\w+)\\s+"([^"]+)', elem)
                        if m:
                            field = m.group(1)
                            value = m.group(2)
                            prop[field] = value
                    #print(prop)
                    trid = prop.get('transcript_id', None)
                    gene_names.add(prop.get('gene_name'))
                    #print(trid, prop)
                                   
                    if trid is None: continue
                    n_tr += 1
                    if trid not in genes:
                        genes[trid] = {'name':prop.get('gene_name', trid),
                                       'chromosome':re.sub('^chr', '', items[0]), 'orientation':items[6], 'exons':[],
                                       'biotype':prop.get('biotype', None)}
                    if items[2] == 'start_codon':
                        genes[trid][['coding_start','coding_end'][int(items[6]=='+')]] = int(items[3])
                    elif items[2] == 'stop_codon':
                        genes[trid][['coding_start','coding_end'][int(items[6]=='-')]] = int(items[3])
                    else:
                        genes[trid]['exons'].append((int(prop.get('exon_number', 0)), int(items[3]), int(items[4])))
                        pass
                    pass
                pass
            #print(n_tr, n_skip, 'n_names', len(gene_names), 'n_genes', len(genes))
            converted = []
            for trid, obj in genes.items():
                g = AnnotatedGene(line=None)
                g.chromosome = obj['chromosome']
                g.strand = obj['orientation']

                g.gene_name = obj['name']
                g.transcript_id = trid
                exons = obj['exons']
                estart = [x_[1] for x_ in exons]
                eend = [x_[2] for x_ in exons]
                g.exon_start = estart
                g.exon_end = eend
                g.start = min(estart)
                g.end = max(eend)
                g.exon_lengths = [(eend[i] - estart[i]) for i in range(len(estart))]
                g.total_length = sum(g.exon_lengths)
                g.coding_start = obj.get('coding_start', -1)
                g.coding_end = obj.get('coding_end', -1)
                g.biotype = obj.get('biotype', None)
                #g.__key = g.__generate_key(g.chromosome, (g.start + g.end) // 2)

                converted.append(g)
            with gzip.open(cache, 'wb') as fo:
                pickle.dump(converted, fo)
                pass
            #print(len(genes), 'transcripts', len(converted), 'genes')
            genes = converted
            pass
        logger.info('{} genes'.format(len(genes)))
        genes = list(sorted(genes, key=functools.cmp_to_key(AnnotatedGene.__comparator)))#lambda x,y:-1 if (x[0]<y[0]) else 1 if (x[0]>y[0]) else -1 if x[1]<y[1] else 1 if x[1]>y[1] else 0)))
    
        return genes
        pass
    
    # @classmethod
    # def __generate_key(cls, chrom, pos):
    #     m = re.match('chr([\\d+XYMTt]+)$', chrom)
    #     if m:
    #         c = m.group(1)
    #         if c.isdigit():
    #             cnum = int(c)
    #         elif c == 'X':
    #             cnum = 100
    #         elif c == 'Y':
    #             cnum = 101
    #         elif c.startswith('M'):
    #             cnum = 200
    #         else:
    #             cnum = 500
    #         return '{:06d}:{:08d}'.format(cnum, pos)
    #     else: # un
    #         return '{}:{:08d}'.format(chrom, pos)
    @classmethod
    def load_genes(cls, filename, **kwargs):
        genes = []
        with open(filename) as fi:
            for line in fi:
                try:
                    gene = AnnotatedGene(line)
                    genes.append(gene)
                except Exception as e:
                    sys.stderr.write(str(e) + ' ' + line + '\n')
                    raise
                    pass
                pass
            pass
        return list(sorted(genes, key=functools.cmp_to_key(__comparator)))
    #lambda x,y:-1 if (x[0]<y[0]) else 1 if (x[0]>y[0]) else -1 if x[1]<y[1] else 1 if x[1]>y[1] else 0)))
    #return list(sorted(genes, key=lambda g:g.position_key))

    @staticmethod
    def __comparator(x, y):
        if x.chromosome != y.chromosome:
            return -1 if x.chromosome < y.chromosome else 1
        if x.position == y.position:
            return 0
        elif x.position < y.position:
            return -1
        return 1

    position = property(lambda x:(x.start + x.end) // 2)
    
    @classmethod
    def find_gene(cls, genes, chromosome, position):
        chromosome = re.sub('^chr', '', chromosome)
        #key = cls.__generate_key(chromosome, position)
        left = 0
        right = len(genes)
        center = (left + right) // 2
        while left < right:
            g = genes[center]
            #print(f'{center}\t{chromosome}:{position}\t{g.chromosome}:{g.start}-{g.end}')
            if g.chromosome < chromosome:
                left = center + 1
            elif g.chromosome > chromosome:
                right = center
                #print(g.__key, left, right, center, chromosome, position)
            elif g.position < position:
                left = center + 1
            elif g.position > position:
                right = center
            else:
                break
            center = (left + right) // 2
            pass
        i = center
        if i < 0 or i >= len(genes):
            return None
        tolerance = 100000
        while i >= 0:
            g = genes[i]
            #print(i, g.__key, g.gene_name, chromosome, position)
            if g.chromosome != chromosome or g.end < position - tolerance:
                break
            else:
                e = g.find_exon(chromosome, position)
                if e:
                    return e
            i -= 1
            pass
        i = center + 1
        while i < len(genes):
            g = genes[i]
            #print(i, g.__key, g.gene_name, chromosome, position)
            if g.chromosome != chromosome or g.start > position + tolerance:
                break
            else:
                e = g.find_exon(chromosome, position)
                if e:
                    return e
                pass
            i += 1
            pass
        return None
    
    position_key = property(lambda s:s.__key)
    pass

def annotate_snps(genes, fn_vcf, fn_out, **kwargs):
    depth = kwargs.get('depth', 1000)
    bias = kwargs.get('bias', 10) * 0.01
    logger = kwargs.get('logger', logging.getLogger())
    allow_xy = kwargs.get('xy', False)
    allow_mt = kwargs.get('mt', False)
    allow_un = kwargs.get('un', False)
    genic_only = kwargs.get('genic', False)
    allow_un = kwargs.get('allow_unassembled', False)
    n_snps = 0
    freq = {}
    gene_symbols = set()
    n_covered = 0
    with open(fn_vcf) as fi, open(fn_out, 'w') as fo:
        for line in fi:
            if line.startswith('#'):
                if not allow_un:
                    m = re.match('##contig=<ID=([^,]+)', line)
                    if m is not None and m.group(1).find('_') >= 0:
                        continue
                fo.write(line)
                continue
            items = line.strip().split('\t')
            if len(items) < 8:
                continue
            if not allow_un and re.match('(chr)?[0-9XYMTm]+$', items[0]) is None:
                continue
            if not allow_xy and (items[0].find('X') >= 0 or items[0].find('Y') >= 0):
                continue
            if not allow_mt and (items[0].find('M') >= 0):
                continue
            m = re.search('\\bDP=(\\d+)', items[7])
            dp = int(m.group(1))
            freq[dp] = freq.get(dp, 0) + 1
            if dp < depth: continue
            n_covered += 1
            m = re.search('\\bDP4=([\\d,]+)', items[7])
            fr, rr, fa, ra = [int(_) for _ in m.group(1).split(',')]
            ratio = (fa + ra) / (fr + rr + fa + ra)
            genic = False
            if bias <= ratio <= 1 - bias:
                chrom = items[0]
                pos = int(items[1])
                mappedgene = AnnotatedGene.find_gene(genes, chrom, pos)
                if mappedgene is not None:
                    gene, enum, bt, offset = mappedgene
                    gene_symbols.add(gene)
                    items[7] = items[7] + ';GENE={}'.format(gene.gene_name)
                    genic = True
            if genic_only is False or genic:
                fo.write('\t'.join(items) + '\n')
                n_snps += 1
            pass
        pass
    return {'n>={}'.format(depth):n_covered, 'n_snps':n_snps, 'n_genes':len(gene_symbols), 'depth_freq':freq, 'n_input':sum(freq.values())}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', default='canfam4/AnnotatedGene.txt')
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('-i', nargs='+')
    parser.add_argument('-o', default='vcf_out')
    parser.add_argument('-d', '--depth', type=int, default=1000)
    parser.add_argument('-b', '--bias', type=int, default=10)
    parser.add_argument('--genic', action='store_true')

    logger = tkutil.get_logger(os.path.basename(__file__))
    args = parser.parse_args()

    verbose = args.verbose
    depth = args.depth
    bias = args.bias
    genic_only = args.genic
    outdir = args.o
    fn_gene = args.g
    fn_info = os.path.join(outdir, 'run.info')
    vcf_files = args.i
    info = {
        'command':sys.argv,
        'gtf':fn_gene, 'genic_only':genic_only, 'bias':bias, 'depth':depth
    }
    if verbose:
        logger.setLevel(10)

    logger.info('loading gene annotation')
    if fn_gene.endswith('.gtf'):
        logger.info('loading GTF file')
        genes = AnnotatedGene.load_gtf(fn_gene, logger=logger)
    else:
        logger.info('loading TSV file')
        genes = AnnotatedGene.load_genes(args.g, logger=logger)
    logger.info('{} genes loaded'.format(len(genes)))
    os.makedirs(outdir, exist_ok=1)
    info['files'] = []
    for vcf in vcf_files:
        fn_out = os.path.join(outdir, os.path.basename(vcf)[:-4] + '.genic.vcf')
        logger.info(f'filtering {vcf} to {fn_out}')
        stat  = annotate_snps(genes, vcf, fn_out, depth=depth, bias=bias, logger=logger, genic=genic_only)
        stat.pop('depth_freq')
        info['files'].append({'input':vcf, 'output':fn_out, 'stat':stat})
        pass
    with open(fn_info, 'w') as fo:
        json.dump(info, fo, indent=2)
                
    
if __name__ == '__main__':
    main()
