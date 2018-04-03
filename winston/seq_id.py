import re

from name_converter import NameConverter

class SeqId:
    """
        Takes dataset contig id
        Extracts internal and external organism names
    """

    SYSTEM_NAME_REGEXP = re.compile(r'(\w{32})_')

    def __init__(self, seqid):
        self.seqid = seqid
        self._prepare_names()

    def is_equal_to(self, other):
        return self.internal_id == other.internal_id and self.external_id == other.external_id

    def _prepare_names(self):
        groups = re.match(SeqId.SYSTEM_NAME_REGEXP, self.seqid).groups()

        if len(groups) == 0:
            raise RuntimeError("Wrong seq_id: %s" % self.seqid)

        self.internal_id = groups[0]
        self.external_id = NameConverter.int_to_ext(self.internal_id)
        self.original_seqid = self.seqid.replace(self.internal_id+'_', '')
