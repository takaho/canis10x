import os, sys, re, argparse
import logging

sys.path.append('/home/data/cellranger/scripts/')

import tkutil
# 590     chr1    689675  702693  TXNL4A.2        100     +       698966  702374  100     5       152,187,70,104,491,     0,519,2607,9231,12527,

class uuGene(object):
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
            return
            
        items = line.strip().split('\t')
        self.chromosome = items[1]
        self.start = int(items[2])
        self.end = int(items[3])
        m = re.match('(.*)\\.(\\d+)$', items[4])
        self.gene_name = m.group(1)
        self.gene_number = int(m.group(2))
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
        self.__key = self.__generate_key(self.chromosome, (self.start + self.end) // 2)
            
    def find_exon(self, chromosome, position):
        n_exons = len(self.exon_start)
        if self.chromosome != chromosome or position < self.start or self.end < position:
            return None
        for i in range(n_exons):
            if self.exon_start[i] <= position <= self.exon_end[i]:
                exon_num = i if self.strand == '+' else n_exons - 1 - i
                is_coding = self.coding_start <= position < self.coding_end
                offset = sum(self.exon_length[:i]) + position - self.exon_start[i]
                if self.strand == '+':
                    offset = self.total_length - offset
                return [self, i, is_coding, offset]
        return None
    def __get_position_key(self):
        return self.__key
    
    @classmethod
    def load_gtf(cls, filename, **kwargs):
        """
GTF
chr1    refGene CDS     15266786        15267037        .       -       2       gene_id "CDH20"; transcript_id "NM_001286988"; exon_number "2"; exon_id "NM_001286988.2"; gene_name "CDH20";
"""
        import gzip, hashlib, pickle
        cache = os.path.join(os.path.dirname(filename), '.' + os.path.basename(filename).split('.')[0])
        forced = kwargs.get('forced', False)
        forced = True
        print(forced)
        if tkutil.check_file(cache, filename) and not forced:
            with gzip.open(cache, 'rb') as fi:
                genes = pickle.loads(fi)
                pass
            pass
        else:
            genes = {}
            n_tr = n_skip = 0
            gene_names = set()
            
            with open(filename) as fi:
                for line in fi:
                    items = line.strip().split('\t')
                    if len(items) < 9 or items[2] != 'exon' or items[3].isdigit() is False:
                        if items[2] == 'exon':
                            print('skip')
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
                        genes[trid] = {'name':prop.get('gene_name', trid), 'chromosome':items[0], 'orientation':items[6], 'exons':[]}
                    genes[trid]['exons'].append((int(prop.get('exon_number', 0)), int(items[3]), int(items[4])))
                    pass
                pass
            print(n_tr, n_skip, 'n_names', len(gene_names), 'n_genes', len(genes))
            converted = []
            for trid, obj in genes.items():
                g = uuGene(line=None)
                g.chromosome = obj['chromosome']
                exons = obj['exons']
                estart = [x_[1] for x_ in exons]
                eend = [x_[2] for x_ in exons]
                g.exon_start = min(estart)
                g.exon_end = max(eend)
                g.exon_lengths = [(eend[i] - estart[i]) for i in range(len(estart))]
                g.total_length = sum(g.exon_lengths)
                g.gene_name = obj['name']
                g.strand = obj['orientation']
                g.__key = g.__generate_key(g.chromosome, (g.start + g.end) // 2)
                converted.append(g)
            with gzip.open(cache, 'wb') as fo:
                pickle.dump(converted, fo)
                pass
            print(len(genes), 'transcripts', len(converted), 'genes')
            genes = converted
            exit()
            pass
        exit()
        return genes
        pass
    
    @classmethod
    def __generate_key(cls, chrom, pos):
        m = re.match('chr([\\d+XYMTt]+)$', chrom)
        if m:
            c = m.group(1)
            if c.isdigit():
                cnum = int(c)
            elif c == 'X':
                cnum = 100
            elif c == 'Y':
                cnum = 101
            elif c.startswith('M'):
                cnum = 200
            else:
                cnum = 500
            return '{:06d}:{:08d}'.format(cnum, pos)
        else: # un
            return '{}:{:08d}'.format(chrom, pos)
    @classmethod
    def load_genes(cls, filename, **kwargs):
        genes = []
        with open(filename) as fi:
            for line in fi:
                try:
                    gene = uuGene(line)
                    genes.append(gene)
                except Exception as e:
                    sys.stderr.write(str(e) + ' ' + line + '\n')
                    raise
                    pass
                pass
            pass
        return list(sorted(genes, key=lambda g:g.position_key))
    @classmethod
    def find_gene(cls, genes, chromosome, position):
        key = cls.__generate_key(chromosome, position)
        left = 0
        right = len(genes)
        center = (left + right) // 2
        while left < right:
            g = genes[center]
            #print(g.__key, left, right, center, chromosome, position)
            if g.__key < key:
                left = center + 1
            elif g.__key > key:
                right = center
            else:
                break
            center = (left + right) // 2
            pass
        i = center
        tolerance = 100000
        while i >= 0:
            g = genes[i]
            print(i, g.__key, g.gene_name, chromosome, position)
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
            print(i, g.__key, g.gene_name, chromosome, position)
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
    n_snps = 0
    with open(fn_vcf) as fi, open(fn_out, 'w') as fo:
        for line in fi:
            if line.startswith('#'):
                fo.write(line)
                continue
            items = line.strip().split('\t')
            #print(items[0])
            if not allow_un and re.match('(chr)?[0-9XYMTm]+$', items[0]) is None:
                continue
            if not allow_xy and (items[0].find('X') >= 0 or items[0].find('Y') >= 0):
                continue
            if not allow_mt and (items[0].find('M') >= 0):
                continue
            m = re.search('\\bDP=(\\d+)', items[7])
            #print(m)
            dp = int(m.group(1))
            if dp < depth: continue
            #print(dp)
            m = re.search('\\bDP4=([\\d,]+)', items[7])
            fr, rr, fa, ra = [int(_) for _ in m.group(1).split(',')]
            ratio = (fa + ra) / (fr + rr + fa + ra)
            genic = False
            #print(fr, rr, fa, ra, ratio)
            if bias <= ratio <= 1 - bias:
                chrom = items[0]
                pos = int(items[1])
                gene = uuGene.find_gene(genes, chrom, pos)
                if gene is not None:
                    items[7] = items[7] + ';GENE={}'.format(gene[0].gene_name)
                    genic = True
                    print('mapped on ', gene)
            if genic_only is False or genic:
                fo.write('\t'.join(items) + '\n')
                n_snps += 1
            pass
        pass
    return n_snps

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', default='canfam4/uuGene.txt')
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
    if verbose:
        logger.setLevel(10)

    logger.info('loading UCSC genes')
    fn_gene = args.g
    if fn_gene.endswith('.gtf'):
        logger.info('loading GTF file')
        genes = uuGene.load_gtf(fn_gene)
    else:
        logger.info('loading TSV file')
        genes = uuGene.load_genes(args.g)
    logger.info('{} genes loaded'.format(len(genes)))
    for vcf in args.i:
        os.makedirs(outdir, exist_ok=1)
        fn_out = os.path.join(outdir, os.path.basename(vcf)[:-4] + '.genic.vcf')
        logger.info(f'filtering {vcf} to {fn_out}')
        annotate_snps(genes, vcf, fn_out, depth=depth, bias=bias, logger=logger, genic=genic_only)
        pass
    
                    
# def load_uugene(filename, **kwargs):
#     verbose = kwargs.get('logger', logging.getLogger())
#     genes = []
#     with open(filename) as fi:
#         for line in fi:
#             items = line.strip().split('\t')
#             chrom = items[1]
#             start = int(items[2])
#             end = int(items[3])
#             m = re.match('(.*)\\.\\d+$', items[4])
#             gene_name = m.group(1)
#             num = int(m.group(2))
#             strand = items[6]
#             coding_start = int(items[7])
#             coding_end = int(items[8])
#             exon_lengths = [int(x_) for x_ in items[11].strip(',').split(',')]
#             exon_start = [int(x_) for x_ in items[11].strip(',').split(',')]
#             #num = gene_name[gene_name.rfind
#                     #genes.append(
                
    
if __name__ == '__main__':
    main()
