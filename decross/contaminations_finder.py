import os, glob

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from decross.blastab.one_vs_all import OneVsAll
from types_manager import TypesManager
from path_resolver import PathResolver
from seq_id import SeqId
from coverage_detector import CoverageDetector
from name_converter import NameConverter
from settings import Settings
from data_manager import Dataset

class ContaminationsFinder:
    """
        Takes one_vs_all blastab file
        Saves clean, contaminated contigs and statistics files
    """

    FILES = {
              'clean': 'fasta',
              'deleted': 'fasta',
              'records': 'csv',
              'deleted_stats': 'csv',
              'missing_kmers': 'csv',
              'contaminations': 'csv' }

    DEFAULT_OWN_KMER = 10**25
    DEFAULT_HIT_KMER = 10**-47

    NO_CONTAMINATION_TYPE = 'NO'

    def __init__(self, file_path):
        self.file_path = file_path
        self.logs = { name: None for name in self.FILES }
        self.external_name = os.path.splitext(os.path.basename(file_path))[0]
        self.internal_name = NameConverter.ext_to_int(self.external_name)
        self.blastab = OneVsAll(file_path)
        self.dataset = Dataset(self.external_name)
        self.contaminations = []

        self.coverage_detector = CoverageDetector()

    def process(self):
        self._open_logs()

        for seq_id, hits in self.blastab.hits_dict.iteritems():
            self._analyze_sequence(SeqId(seq_id), hits)

        self._log_contaminations()
        self._log_results()

        self._close_logs()

    def _analyze_sequence(self, seq_id, hits):
        own_kmer = self._detect_kmer(seq_id.seqid, self.DEFAULT_OWN_KMER)

        for hit in hits:
            hit_kmer = self._detect_kmer(hit.subject_seq_id.seqid, self.DEFAULT_HIT_KMER)
            self.logs['records'].write('%s,%s,%s,%s\n' % (seq_id.original_seqid,
                                                          hit.subject_seq_id.original_seqid,
                                                          own_kmer,
                                                          hit_kmer))

            # check type and threshold
            pair_type = TypesManager.get_type(seq_id.external_id, hit.subject_seq_id.external_id)

            if pair_type == self.NO_CONTAMINATION_TYPE: # NO
                continue

            ratio = getattr(Settings.decross.coverage_ratio, pair_type.lower())

            if not ratio:
                raise Exception("Cannot detect ratio")

            if hit_kmer == 0:
                hit_kmer = self.DEFAULT_HIT_KMER

            if own_kmer/hit_kmer <= ratio: # our is less or equal than 1.5x of their
                self.contaminations.append(hit)

                threshold = TypesManager.get_threshold(hit.query_seq_id.external_id,
                                                       hit.subject_seq_id.external_id)

                line = (hit.query_seq_id.original_seqid, hit.subject_seq_id.original_seqid, hit.pident,
                        hit.calculate_qlen(), hit.qcovhsp, own_kmer,
                        hit_kmer, pair_type, threshold)
                self.logs['deleted_stats'].write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % line)

    def _open_logs(self):
        for name, ext in self.FILES.iteritems():
            f_name = '%s_%s.%s' % (self.external_name, name, ext)
            self.logs[name] = open(PathResolver.results_path_for(f_name), 'w')

    def _close_logs(self):
        for name in self.logs:
            self.logs[name].close()

    def _dataset_path(self):
        return self.dataset.contigs_output_path()

    def _detect_kmer(self, contig_id, default_value):
        coverage = self.coverage_detector.kmer_by_contig_id(contig_id)

        if coverage != None:
            return float(coverage)
        else:
            self.logs['missing_kmers'].write('%s\n' % (contig_id))
            return float(default_value)

    def _log_contaminations(self):
        contaminated_from = [h.subject_seq_id.external_id for h in self.contaminations]
        stats = { i: contaminated_from.count(i) for i in list(set(contaminated_from)) }

        for name, cnt in stats.iteritems():
            self.logs['contaminations'].write('%s,%s\n' % (name, cnt))

    def _log_results(self):
        contaminated_ids = [s.query_seq_id.seqid for s in self.contaminations]

        for record in SeqIO.parse(self._dataset_path(), "fasta"):
            current_seq_id = record.id
            record.id = record.description = SeqId(record.id).original_seqid

            if current_seq_id in contaminated_ids:
                SeqIO.write(record, self.logs['deleted'], "fasta")
            else:
                SeqIO.write(record, self.logs['clean'], "fasta")
