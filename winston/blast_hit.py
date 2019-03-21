from seq_id import SeqId
from coverage_detector import CoverageDetector

class BlastHit():
    LEN_RATIO = 70
    LEN_MINIMUM = 100

    DEFAULT_QUERY_RPKM = 10 ** 25
    DEFAULT_SUBJECT_RPKM = 10 ** -47

    def __init__(self, raw_hit):
        self.raw = raw_hit
        splitted = raw_hit.strip().split()

        self.qseqid = splitted[0]
        self.sseqid = splitted[1]
        self.pident = float(splitted[2])
        self.length = int(splitted[3])
        self.qstart = int(splitted[7])
        self.qend = int(splitted[8])
        self.qcovhsp = int(splitted[-1])

        self.query_seq_id = SeqId(self.qseqid)
        self.subject_seq_id = SeqId(self.sseqid)

        self._query_RPKM = None
        self._subject_RPKM = None
        self.missing_RPKMs = []

    def calculate_qlen(self):
        return self.qend - (self.qstart - 1)

    def is_reliable(self):
        return self.qcovhsp >= self.LEN_RATIO and self.calculate_qlen() >= self.LEN_MINIMUM

    def equal_organisms(self):
        return self.query_seq_id.is_equal_to(self.subject_seq_id)

    def query_RPKM(self):
        if not self._query_RPKM:
            self._query_RPKM = self._detect_RPKM(self.query_seq_id.seqid)

            if self._query_RPKM == None:
                self.missing_RPKMs.append(self.query_seq_id.original_seqid)
                self._query_RPKM = self.DEFAULT_QUERY_RPKM

        return self._query_RPKM

    def subject_RPKM(self):
        if not self._subject_RPKM:
            self._subject_RPKM = self._detect_RPKM(self.subject_seq_id.seqid)

            if self._subject_RPKM == None:
                self.missing_RPKMs.append(self.subject_seq_id.original_seqid)
                self._subject_RPKM = self.DEFAULT_SUBJECT_RPKM

        return self._subject_RPKM

    def _detect_RPKM(self, contig_id):
        rpkm = CoverageDetector().kmer_by_contig_id(contig_id)

        if rpkm != None:
            return float(rpkm)
        else:
            return None

    def rpkm_difference(self):
        return abs(self.query_RPKM() - self.subject_RPKM())

    def to_s(self):
        line = (self.query_seq_id.original_seqid, self.subject_seq_id.original_seqid, self.pident,
                self.calculate_qlen(), self.qcovhsp, self.query_RPKM(), self.subject_RPKM())
        return "%s,%s,%s,%s,%s,%s,%s" % line

    @staticmethod
    def generate(file_path):
        return [BlastHit(raw_hit) for raw_hit in open(file_path).readlines()]
