import os, glob

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from winston.blastab.one_vs_all import OneVsAll
from winston.blast_hit import BlastHit
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
              'suspicious_hits': 'csv',
              'missing_coverage': 'csv',
              'contaminations': 'csv',
              'contamination_sources': 'csv'}

    NO_CONTAMINATION_TYPE = 'NO'

    def __init__(self, file_path):
        self.file_path = file_path
        self.logs = { name: None for name in self.FILES }
        self.external_name = os.path.splitext(os.path.basename(file_path))[0]
        self.internal_name = NameConverter.ext_to_int(self.external_name)
        self.blastab = OneVsAll(file_path)
        self.dataset = Dataset(self.external_name)

        self.suspicious_hits = {}
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
        own_seq_id = seq_id.seqid
        query_rpkm = hits[0].query_RPKM()

        for hit in hits:
            hit._query_RPKM = query_rpkm # caching
            subject_rpkm = hit.subject_RPKM()

            if len(hit.missing_RPKMs) > 0:
                for seqid in hit.missing_RPKMs:
                    self.logs['missing_coverage'].write("%s\n" % (seqid))

            # check type and threshold
            pair_type = TypesManager.get_type(seq_id.external_id, hit.subject_seq_id.external_id)

            if pair_type == self.NO_CONTAMINATION_TYPE: # NO
                continue

            ratio = getattr(Settings.winston.coverage_ratio, pair_type)

            if not ratio:
                raise Exception("Cannot detect ratio")

            if subject_rpkm == 0:
                subject_rpkm = BlastHit.DEFAULT_SUBJECT_RPKM

            # our is less or equal than 1.5x of their
            if query_rpkm/subject_rpkm <= ratio:
                if own_seq_id not in self.suspicious_hits:
                    self.suspicious_hits[own_seq_id] = []

                self.suspicious_hits[own_seq_id].append(hit)

                threshold = TypesManager.get_threshold(hit.query_seq_id.external_id,
                                                       hit.subject_seq_id.external_id)

                hit_line = hit.to_s()
                self.logs['suspicious_hits'].write("%s,%s,%s\n" % (hit_line, pair_type, threshold))

        if own_seq_id in self.suspicious_hits and len(self.suspicious_hits[own_seq_id]) > 0:
            self.contaminations.append(self._get_best_hit(self.suspicious_hits[own_seq_id]))

    def _open_logs(self):
        for name, ext in self.FILES.iteritems():
            f_name = '%s_%s.%s' % (self.external_name, name, ext)
            self.logs[name] = open(PathResolver.results_path_for(f_name), 'w')

    def _close_logs(self):
        for name in self.logs:
            self.logs[name].close()

    def _dataset_path(self):
        return self.dataset.contigs_output_path()

    def _log_contaminations(self):
        sources = {}

        for hit in self.contaminations:
            self.logs['contaminations'].write('%s\n' % hit.to_s())

            from_id = hit.subject_seq_id.external_id

            if from_id not in sources:
                sources[from_id] = 1
            else:
                sources[from_id] += 1

        for from_id, cnt in sources.items():
            self.logs['contamination_sources'].write('%s,%s\n' % (from_id, cnt))


    def _log_results(self):
        contaminated_ids = [s.query_seq_id.seqid for s in self.contaminations]

        for record in SeqIO.parse(self._dataset_path(), "fasta"):
            current_seq_id = record.id
            record.id = record.description = SeqId(record.id).original_seqid

            if current_seq_id in contaminated_ids:
                SeqIO.write(record, self.logs['deleted'], "fasta")
            else:
                SeqIO.write(record, self.logs['clean'], "fasta")

    def _get_best_hit(self, hits):
        best_hit = None

        for hit in hits:
            if not best_hit:
                best_hit = hit
                continue
            if hit.rpkm_difference() > best_hit.rpkm_difference():
                best_hit = hit

        return best_hit
