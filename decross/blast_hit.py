from seq_id import SeqId

class BlastHit():
    LEN_RATIO = 70
    LEN_MINIMUM = 100

    def __init__(self, raw_hit):
        splitted = raw_hit.strip().split()

        if len(splitted) < 3:
            print splitted

        self.qseqid = splitted[0]
        self.sseqid = splitted[1]
        self.pident = float(splitted[2])
        self.length = int(splitted[3])
        self.qstart = int(splitted[7])
        self.qend = int(splitted[8])
        self.qcovhsp = int(splitted[-1])

        self.query_seq_id = SeqId(self.qseqid)
        self.subject_seq_id = SeqId(self.sseqid)

    def calculate_qlen(self):
        return self.qend - (self.qstart - 1)

    def is_reliable(self):
        return self.qcovhsp >= self.LEN_RATIO and self.calculate_qlen() >= self.LEN_MINIMUM

    def equal_organisms(self):
        return self.query_seq_id.is_equal_to(self.subject_seq_id)

    @staticmethod
    def generate(file_path):
        return [BlastHit(raw_hit) for raw_hit in open(file_path).readlines()]
