#!/usr/bin/env python3

import glob, os, re
import ntpath
import argparse
import subprocess
from tqdm import tqdm


description = 'Get global statistics of contaminations for the Winston output'
parser = argparse.ArgumentParser(description=description)
parser.add_argument("-r", "--results_path", required=True, help="Winston results folder")
parser.add_argument("-o", "--output", required=True, help="Full statistics output file (.csv)")
parser.add_argument("-t", "--types", required=False, help="Comma separated types for calculation, all if not provided")
options = parser.parse_args()

def get_contigs_count(path):
    result = 0

    if ntpath.exists(path):
        try:
            result = int(subprocess.check_output(f"grep -c '>' {path}", shell=True))
        except subprocess.CalledProcessError as grepexc:
            pass

    return result


def main():
    paths = glob.glob(ntpath.join(options.results_path, "*_deleted_stats.csv"))
    paths.sort()
    acc_regexp = re.compile("MMETSP\d{4}")

    types = []
    if options.types:
        types = options.types.split(',')

    without_type_restriction = len(types) == 0

    with open(options.output, 'w') as out_f:
        out_f.write("acc_id,contig_cnt,deleted_cnt\n")

        for path in tqdm(paths):
            acc_id = acc_regexp.search(ntpath.basename(path))[0]
            deleted_cnt = 0

            with open(path) as f:
                last_contig = None

                for line in f.readlines():
                    data = line.strip().split(',')
                    contig = data[0]

                    if without_type_restriction or data[-2] in types:
                        if contig != last_contig:
                            deleted_cnt += 1
                            last_contig = contig

            clean_path = ntpath.join(options.results_path, f"{acc_id}_clean.fasta")
            deleted_path = ntpath.join(options.results_path, f"{acc_id}_deleted.fasta")

            all_clean_cnt = get_contigs_count(clean_path)
            all_deleted_cnt = get_contigs_count(deleted_path)

            out_f.write(f"{acc_id},{all_clean_cnt+all_deleted_cnt},{deleted_cnt}\n")


if __name__ == '__main__':
    main()
