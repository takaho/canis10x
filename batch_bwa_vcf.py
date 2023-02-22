import os, sys, re
import subprocess
import argparse

def check_file_is_processed(filename, date_stat=None):
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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--genome', default='../cellranger_db/BWA/canFam4.fa', help='Genome fasta')
    parser.add_argument('-d', '--db', default='../cellranger_db/BWA/canFam4.fa', help='BWA database')
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--test', action='store_true')
    parser.add_argument('-p', default=4, type=int, help='number of threads')
    parser.add_argument('-i', nargs='+', default='fastq')
    parser.add_argument('-o', default='out')
    parser.add_argument('--paired', action='store_true')

    args = parser.parse_args()
    outdir = args.o
    filenames = args.i
    os.makedirs(outdir, exist_ok=1)

    test = args.test
    n_threads = args.p
    verbose = args.verbose
    paired = args.paired
    bwadb = args.db#'/home/data/endo/canis/cellranger_db/ROSCfam1/fasta/ROSCfam1'
    fn_genome = args.genome# '/home/data/endo/canis/cellranger_db/ROSCfam1/fasta/genome.fa'

    #cmd = 'bwa', 'mem', '-t', '16', '-o', fn_sam, bwadb, fn_fasta

    dataset = {}


    if paired:
        pat = re.compile('(\\w+)_R(1|2)(_\\d{3})?\\.fastq.gz$')
    else:
        pat = re.compile('(\\w+)?_.*\\.fastq.gz$')

    for fn in filenames:
        if fn.endswith('.sam'):
            name = os.path.basename(fn).split('.')[0]
            fn_bam = os.path.join(outdir, f'{name}.bwa.bam')
            fn_vcf = os.path.join(outdir, f'{name}.bwa.vcf')
            dataset[name] = [None, None, fn_sam, fn_bam, fn_vcf]
        elif fn.endswith('.bam'):
            name = os.path.basename(fn).split('.')[0]
            fn_bam = fn
            fn_vcf = os.path.join(outdir, f'{name}.bwa.vcf')
            dataset[name] = [None, None, None, fn_bam, fn_vcf]
        else:
            m = pat.match(os.path.basename(fn))#re.match('(\\w+)_R(1|2)(_\\d{3})?\\.fastq.gz$', fn)
            if m is None:
               if verbose:
                   sys.stderr.write(f'skip {fn}\n')
                   continue
            name = m.group(1)
            if name not in dataset:
                fn_sam = os.path.join(outdir, f'{name}.bwa.sam')
                fn_bam = os.path.join(outdir, f'{name}.bwa.bam')
                fn_vcf = os.path.join(outdir, f'{name}.bwa.vcf')
                dataset[name] = [None, None, fn_sam, fn_bam, fn_vcf]
            if paired:
                r = int(m.group(2))
                if r == 1:
                    dataset[name][0] = fn
                elif r == 2:
                    dataset[name][1] = fn
                    pass
                pass
            else:
                dataset[name][0] = fn

    for name in dataset.keys():
        r1, r2, fn_sam, fn_bam, fn_vcf = dataset[name]
        if paired:
            if r1 is not None or r2 is not None or sam is None:
                if verbose:
                    sys.stderr.write('sample {} is incomplet\n')
                if not check_file_is_processed(fn_sam) and not check_file_is_processed(fn_bam):
                    continue
                pass
            pass
        # mapping
        if r1 is not None and check_file_is_processed(fn_bam, r1) is False and check_file_is_processed(fn_sam, r1) is False:
            if paired:
                cmd = 'bwa', 'mem', '-t', str(n_threads), '-o', fn_sam, bwadb, r1, r2
            else:
                cmd = 'bwa', 'mem', '-t', str(n_threads), '-o', fn_sam, bwadb, r1 or r2
            if verbose:
                sys.stderr.write(' '.join(cmd) + '\n')
            if not test:
                proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
                with open(fn_sam, 'wb') as fo:
                    for line in proc.stdout:
                        fo.write(line)
                proc.stdout.close()
                proc.wait()
        else:
            if verbose:
                sys.stderr.write(f'{fn_sam} is already processed\n')
                pass
            pass
        
        if fn_sam is not None and check_file_is_processed(fn_bam, fn_sam) is False:
            cmd = 'samtools', 'sort', '-@', '16', '-o', fn_bam, fn_sam
            if verbose:
                sys.stderr.write(' '.join(cmd) + '\n')
            if not test:
                subprocess.Popen(cmd).wait()
                if check_file_is_processed(fn_bam, fn_sam):
                    os.unlink(fn_sam)
                    pass
                pass
            pass
        elif verbose:
            sys.stderr.write(f'{fn_bam} is already processed\n')
            pass
        fn_bai = fn_bam + '.bai'
        if check_file_is_processed(fn_bai, fn_bam) is False:
            cmd = 'samtools', 'index', '-@', '16', fn_bam
            if verbose: sys.stderr.write(' '.join(cmd) + '\n')
            if not test:
                subprocess.Popen(cmd).wait()
                pass
            pass
        elif verbose:
            sys.stderr.write('already indexed\n')
            pass
        if check_file_is_processed(fn_vcf, fn_bam) is False:
            cmd = 'samtools', 'mpileup', '-Ou', '-f', fn_genome, fn_bam, '|', 'bcftools', 'call', '-mv', '-Ob', '-o', fn_vcf
            #../cellranger_db/BWA/canFam4.fa test/S1_S1_R2.bwa.bam | bcftools call -mv -Ob -o test/S1_S1_R2.vcf/S1_S1_R2.bwa.vcf
            #cmd = 'samtools', 'mpileup', '-uf', fn_genome, fn_bam, '|', 'bcftools', 'view', '-mv', '-o', fn_vcf#'-', '>', fn_vcf# hs37d5_allseqs_bwa.raw.vcf
            if verbose:
                sys.stderr.write(' '.join(cmd) + '\n')
            if not test:
                subprocess.Popen(cmd).wait()
                pass
            pass
        elif verbose:
            sys.stderr.write(f'vcf is alread generated.')
            pass
        pass


if __name__ == '__main__':
    main()
