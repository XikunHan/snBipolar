# -*- coding: utf-8 -*-
###############################################################################
# snBD — PLINK clumping of GWAS or eQTL summary statistics
#
# Purpose:
#   Identify LD-clumped variants and export their summary-statistic rows.
#
# Inputs:
#   Required: --dir_path, --gwas_summary_file, --trait.
#   Optional SNP/P column names, PLINK LD prefix and clumping thresholds.
#   Set the PLINK environment variable to specify the executable.
#
# Outputs:
#   PLINK <trait>_GWAS_result_all.* files and an indep_GWhits_<trait>_... file
#   whose name records the window, r2, p1 and p2 thresholds.
#
# Method overview:
#   Run PLINK clumping, extract the lead-variant rows and sort the output.
#   Default settings: 1000 kb, r2=0.01, p1=p2=5e-8; reference MAF >=0.01.
#
# Dependencies:
#   Python standard library; PLINK1.9, awk, sed, sort, zcat.
#
# Usage:
#   From repository root: python mendelian_randomization/clump_gwas.py --help
#
###############################################################################

import argparse
import os
import subprocess
import gzip
import shlex
import tempfile


# Function: is_gz_file()
# Purpose: identify gzip input from its two-byte signature, not its extension.
# Inputs: filepath (existing input file).
# Returns: True for a gzip signature; does not decompress the file.
def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'
        
# Function: clump_gwas()
# Purpose: run PLINK clumping and select matching summary-statistic rows.
# Inputs: parsed command-line arguments listed at the bottom of this script.
# Outputs: PLINK files and indep_GWhits_<trait>_...; returns None.
# Shell arguments are quoted; a failed command stops the pipeline.
def clump_gwas(args):
    if args.dir_path is None:
        raise ValueError('The --dir_path flag is required.')
    if args.gwas_summary_file is None:
        raise ValueError('The --gwas_summary_file flag is required.')
    if args.trait is None:
        raise ValueError('The --trait flag is required.')
    
    
    dir_path = args.dir_path
    os.makedirs(dir_path, exist_ok=True)
    gwas_summary_file = args.gwas_summary_file
    extract_file=args.extract
    
    gwas_summary_file_full_path = dir_path + '/'+ gwas_summary_file if os.path.basename(gwas_summary_file) == gwas_summary_file else gwas_summary_file
    out_file = args.dir_path +'/'+ args.trait + '_GWAS_result_all'
    
    independent_file = args.dir_path + '/indep_GWhits_' + args.trait + '_' + args.kbd + 'kb.r2_' + args.r2level + '_p1_' + args.p1thres + '_p2_' + args.p2thres
    for previous in (out_file + '.clumped', independent_file):
        if os.path.exists(previous):
            os.remove(previous)

    plink = os.environ.get("PLINK", "./software/plink_linux_x86_64_20201019_v1.9/plink")

    command = [plink, '--bfile', args.LDref_bfile, '--maf', '0.01',
               '--clump', gwas_summary_file_full_path, '--clump-snp-field', args.clump_snp,
               '--clump-field', args.clump_field, '--clump-p1', args.p1thres,
               '--clump-p2', args.p2thres, '--clump-kb', args.kbd, '--clump-r2', args.r2level,
               '--out', out_file, '--threads', args.threads, '--memory', args.memory]
    if extract_file != 'ALL':
        command += ['--extract', extract_file]
    print(' '.join(shlex.quote(arg) for arg in command))
    subprocess.run(command, check=True)

    clump_OK = os.path.exists(out_file+'.clumped')
    clumped_file=out_file+'.clumped'
    
    
    gz = is_gz_file(gwas_summary_file_full_path)
    
    if gz:
        with gzip.open(gwas_summary_file_full_path, 'r') as f:
            try:
                f_header = f.readline().decode("utf-8").split()
                print('input_file is gzipped')
            except gzip.BadGzipFile:
                f_header = f.readline().split()
    else:
        with open(gwas_summary_file_full_path, 'r') as f:
            f_header = f.readline().split()
    
    print(f_header)
    
    snp_index = f_header.index(args.clump_snp) + 1
    print('Column ' + str(snp_index) + ' is ' + args.clump_snp)
    
    
    if not clump_OK:
        print('PLINK created no .clumped file; inspect its log for the no-hit reason.')
    else:
        print('Clump successed.')
        tmp_file = None
        try:
            selection_input = gwas_summary_file_full_path
            if gz:
                fd, tmp_file = tempfile.mkstemp(prefix='clump_', suffix='.txt', dir=dir_path)
                os.close(fd)
                bashcmd = 'zcat ' + shlex.quote(gwas_summary_file_full_path) + ' > ' + shlex.quote(tmp_file)
                subprocess.run(['bash', '-o', 'pipefail', '-c', bashcmd], check=True)
                selection_input = tmp_file
            awk_program = "NR==FNR{a[$3]; next;}{if(FNR == 1 || $" + str(snp_index) + " in a) print $0}"
            bashcmd = ('awk ' + shlex.quote(awk_program) + ' ' + shlex.quote(clumped_file) +
                       ' ' + shlex.quote(selection_input) + ' | (sed -u 1q; sort -V) > ' +
                       shlex.quote(independent_file))
            print(bashcmd)
            subprocess.run(['bash', '-o', 'pipefail', '-c', bashcmd], check=True)
        except Exception:
            if os.path.exists(independent_file):
                os.remove(independent_file)
            raise
        finally:
            if tmp_file is not None and os.path.exists(tmp_file):
                os.remove(tmp_file)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Clump GWAS summary statistics.')
    parser.add_argument('--dir_path', type=str, required=True, help='path to the folder for GWAS summary statistics.')
    parser.add_argument('--gwas_summary_file', type=str, required=True, help='file name of GWAS summary statistics.')
    parser.add_argument('--trait', type=str, required=True, help='trait name')
    parser.add_argument('--clump_snp', type=str, default='SNP', help='snp column name in GWAS summary statistics')
    parser.add_argument('--clump_field', type=str, default='P', help='P value column name in GWAS summary statistics')
    parser.add_argument('--LDref_bfile', type=str, default='./data/1kg_LD_ref/g1000_eur', help='LD referebce file')
    parser.add_argument('--kbd', type=str, default='1000', help='clump parameter, kb')
    parser.add_argument('--r2level', type=str, default='0.01', help='clump parameter, r2')
    parser.add_argument('--p1thres', type=str, default='5e-8', help='clump parameter, p1 thres')
    parser.add_argument('--p2thres', type=str, default='5e-8', help='clump parameter, p2 thres')
    parser.add_argument('--extract', type=str, default='ALL', help='clump parameter, p2 thres')
    parser.add_argument('--threads', type=str, default='1', help='threads in plink (v1.9)')
    parser.add_argument('--memory', type=str, default='4000', help='memory in plkink (v1.9)')
    
    args = parser.parse_args()
    
    clump_gwas(args)
    
