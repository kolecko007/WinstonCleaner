#!/usr/bin/env python2.7

import glob, os, re
from pathlib import Path
from optparse import OptionParser
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = OptionParser(description="Get statistics for the dataset")
parser.add_option("--output_path", help="Decross output folder to save the results (folder)")
parser.add_option("--dirty_path", help="List of dirty sequences (.txt)")
parser.add_option("--org_id", help="Organism id")
parser.add_option("--for_stats", action="store_true", help="Statistics for table file")
(options, args) = parser.parse_args()

dataset_path = os.path.join(options.output_path,
                           'datasets',
                           "%s.fasta" % options.org_id)
deleted_path = os.path.join(options.output_path,
                           'results',
                           "%s_deleted.fasta" % options.org_id)
clean_path = os.path.join(options.output_path,
                          'results',
                          "%s_clean.fasta" % options.org_id)

def read_from_fasta(path):
  result = []
  for line in open(path, 'r').readlines():
    if line[0] == '>':
      result.append(line.split(' ')[0][1:].strip())
  return result

def read_from_list(path):
  result = []
  for line in open(path, 'r').readlines():
    result.append(line.strip())
  return result

seq_count = int(subprocess.check_output(['grep', '-c', '>', dataset_path]).strip())

dirty = read_from_fasta(dataset_path)
deleted = read_from_fasta(deleted_path)
clean = read_from_fasta(clean_path)
to_delete = read_from_list(options.dirty_path)

should_be_deleted_cnt = len(to_delete)
should_be_kept_cnt = seq_count-len(to_delete)

print options.org_id
print
print "Seq count ", seq_count
print
print "Should be deleted: ", should_be_deleted_cnt
print "Should be kept: ", should_be_kept_cnt
print
print "Actually deleted: ", len(deleted)
print "Actually kept: ", len(clean)
print

correct = list(set(deleted).intersection(set(to_delete)))
print "Correct deletions: ", len(correct)

mistakenly_kept = list(set(to_delete) - set(deleted))
mistakenly_kept_pct = round((len(mistakenly_kept)/float(should_be_deleted_cnt))*100, 2)
print "Mistakenly kept: %s (%s%%)" % (len(mistakenly_kept), mistakenly_kept_pct)

mistakenly_deleted = list(set(deleted) - set(to_delete))
mistakenly_deleted_pct = round((len(mistakenly_deleted)/float(should_be_kept_cnt))*100, 2)
print "Mistakenly deleted: %s (%s%%)" % (len(mistakenly_deleted), mistakenly_deleted_pct)

## Saving statistics to file

mis_deleted_path = os.path.join(options.output_path,
                                'results',
                                "%s_mistakenly_deleted.csv" % options.org_id)

deleted_stats_path = os.path.join(options.output_path,
                                  'results',
                                  "%s_deleted_stats.csv" % options.org_id)

with open(mis_deleted_path, 'w') as f:
    for line in open(deleted_stats_path, 'r').readlines():
      org = line.split(',')[0]
      if org in mistakenly_deleted:
        f.write(line)

## make mis. deleted fasta

mis_deleted_fasta_path = os.path.join(options.output_path,
                                      'results',
                                      "%s_mistakenly_deleted.fasta" % options.org_id)

deleted_fasta_path = os.path.join(options.output_path,
                                  'results',
                                  "%s_deleted.fasta" % options.org_id)

with open(mis_deleted_fasta_path, 'w') as out_f:
  for record in SeqIO.parse(deleted_fasta_path, "fasta"):
    if record.id in mistakenly_deleted:
      SeqIO.write(record, out_f, "fasta")



mis_kept_path = os.path.join(options.output_path,
                            'results',
                            "%s_mistakenly_kept.csv" % options.org_id)

kept_stats_path = os.path.join(options.output_path,
                               'results',
                               "%s_kept_stats.csv" % options.org_id)

all_kept = {}
for line in open(kept_stats_path, 'r').readlines():
  s = line.split(',')
  all_kept[s[0]] = line

with open(mis_kept_path, 'w') as f:
  for kept_contig in mistakenly_kept:
    if kept_contig in all_kept:
      f.write(all_kept[kept_contig])
    else:
      f.write(kept_contig + '\n')


if options.for_stats:
  print
  print options.org_id
  print
  print seq_count
  print
  print should_be_deleted_cnt
  print should_be_kept_cnt
  print
  print len(deleted)
  print len(clean)
  print
  print len(correct)
  print len(mistakenly_kept)
  print len(mistakenly_deleted)
  print mistakenly_kept_pct
  print mistakenly_deleted_pct
