from ..blast_hit import BlastHit
from ..types_manager import TypesManager
from ..settings import Settings
from blastab import Blastab

class OneVsAll(Blastab):
    def __init__(self, file_path):
        super(OneVsAll, self).__init__(file_path)
        self.hits_dict = self._make_hits_dict()

    def _make_hits_dict(self):
        # { left_seq: [<BlastHit object>, ...], ... }
        pair_dict = self._build_pair_dict()

        result = {}
        len_ratio = Settings.winston.hits_filtering.len_ratio
        len_minimum = Settings.winston.hits_filtering.len_minimum

        for pair, blast_hit in pair_dict.iteritems():
            left_seq_id = blast_hit.query_seq_id
            right_seq_id = blast_hit.subject_seq_id

            threshold = TypesManager.get_threshold(left_seq_id.external_id, right_seq_id.external_id)

            if not threshold:
                raise RuntimeError("threshold for %s - %s not found " % (left_seq_id.external_id, right_seq_id.external_id))

            if float(blast_hit.pident) >= threshold:
                if int(blast_hit.qcovhsp) >= len_ratio or int(blast_hit.calculate_qlen()) >= len_minimum:
                    if left_seq_id.seqid not in result:
                        result[left_seq_id.seqid] = [blast_hit]
                    else:
                        result[left_seq_id.seqid].append(blast_hit)

        return result

    def _build_pair_dict(self):
        # { (left_seq, right_seq): <BlastHit object>, ... }
        pair_dict = {}

        for hit in self.hits:
            key = (hit.qseqid, hit.sseqid)
            if not hit.equal_organisms() and key not in pair_dict:
                pair_dict[key] = hit

        return pair_dict
