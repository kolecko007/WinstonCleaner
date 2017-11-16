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

    FILES = { 'deleted': 'fasta',
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

        self.seq_dict = self._make_seq_dict()
        self.contaminations = []
        self.coverage_detector = CoverageDetector()

    def process(self):
        self._open_logs()

        for seq_id, hits in self.blastab.hits_dict.iteritems():
            seq_id = SeqId(seq_id)
            status = self._analyze_sequence(seq_id, hits)

            if status == 'good':
                pass
            elif status == 'bad':
                seq_rec = SeqRecord(Seq(self.seq_dict[seq_id.seqid]), id=seq_id.original_seqid,
                                                                      description='')
                SeqIO.write(seq_rec, self.logs['deleted'], "fasta")
            else:
                raise Exception("Wrong status: %s" % status)

        contaminated_from = [h.subject_seq_id.external_id for h in self.contaminations]
        stats = { i: contaminated_from.count(i) for i in list(set(contaminated_from)) }

        for name, cnt in stats.iteritems():
            self.logs['contaminations'].write('%s,%s\n' % (name, cnt))

        self._close_logs()
        self._clean_dataset()

    def _analyze_sequence(self, seq_id, hits):
        status = 'good'

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
                status = 'bad'

                self.contaminations.append(hit)

                threshold = TypesManager.get_threshold(hit.query_seq_id.external_id,
                                                       hit.subject_seq_id.external_id)

                line = (hit.query_seq_id.original_seqid, hit.subject_seq_id.original_seqid, hit.pident,
                        hit.calculate_qlen(), hit.qcovhsp, own_kmer,
                        hit_kmer, pair_type, threshold)
                self.logs['deleted_stats'].write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % line)

        return status

    def _open_logs(self):
        for name, ext in self.FILES.iteritems():
            f_name = '%s_%s.%s' % (self.external_name, name, ext)
            self.logs[name] = open(PathResolver.results_path_for(f_name), 'w')

    def _close_logs(self):
        for name in self.logs:
            self.logs[name].close()

    def _make_seq_dict(self):
        with open(self._dataset_path(), 'r') as f:
            line = f.read()

        seqs = line.split('>')
        seqs = seqs[1:]
        d = {}
        for seq in seqs:
            d[seq.split()[0].strip()] = ''.join(seq.split('\n')[1:])
        return d

    def _dataset_path(self):
        return self.dataset.contigs_output_path()

    def _detect_kmer(self, contig_id, default_value):
        coverage = self.coverage_detector.kmer_by_contig_id(contig_id)

        if coverage != None:
            return float(coverage)
        else:
            self.logs['missing_kmers'].write('%s\n' % (contig_id))
            return float(default_value)

    def _clean_dataset(self):
        if len(self.contaminations) == 0:
            return True

        in_path = self._dataset_path()
        out_path = PathResolver.results_path_for("%s_clean.fasta" % self.external_name)
        contaminated_ids = [s.query_seq_id.seqid for s in self.contaminations]

        with open(out_path, 'w') as out_f:
            for record in SeqIO.parse(in_path, "fasta"):
                if record.id not in contaminated_ids:
                    record.id = record.description = SeqId(record.id).original_seqid
                    SeqIO.write(record, out_f, "fasta")
