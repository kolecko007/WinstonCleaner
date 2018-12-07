#!/usr/bin/env python3

import glob, os, re
import ntpath
import argparse

import subprocess

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

description = 'Get global statistics of contaminations for the Winston output'
parser = argparse.ArgumentParser(description=description)
parser.add_argument("-r", "--results_path", required=True, help="Winston results folder")
parser.add_argument("-o", "--output", required=True, help="Full statistics output file (.csv)")
parser.add_argument("-t", "--types_path", required=True, help="Path to types file (.csv)")
options = parser.parse_args()


def get_contigs_count(path):
    return int(subprocess.check_output(f"grep -c '>' {path}", shell=True))

def parse_types(path):
    result = {}
    with open(path) as f:
        for line in f.readlines():
            splitted = line.strip().split(',')
            result[(splitted[0], splitted[1])] = splitted[-1]
    return result

def main():
    types = parse_types(options.types_path)
    paths = glob.glob(ntpath.join(options.results_path, "*_contaminations.csv"))
    paths.sort()
    acc_regexp = re.compile("MMETSP\d{4}")

    with open(options.output, 'w') as out_f:
        out_f.write("to_id,from_id,contam_cnt,contig_cnt,clean_cnt,deleted_cnt,type\n")

        for path in paths:
            acc_id = acc_regexp.search(ntpath.basename(path))[0]
            print(acc_id)

            clean_path = ntpath.join(options.results_path, f"{acc_id}_clean.fasta")
            deleted_path = ntpath.join(options.results_path, f"{acc_id}_deleted.fasta")

            clean_cnt = get_contigs_count(clean_path)
            deleted_cnt = get_contigs_count(deleted_path)

            with open(path) as f:
                for line in f.readlines():
                    right_acc_id = line.strip().split(',')[0]
                    out_f.write(f"{acc_id},{line.strip()},{clean_cnt + deleted_cnt},{clean_cnt},{deleted_cnt},{types[(acc_id, right_acc_id)]}\n")

if __name__ == '__main__':
    main()
