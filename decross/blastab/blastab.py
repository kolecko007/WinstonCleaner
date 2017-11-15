from ..blast_hit import BlastHit

class Blastab(object):
    def __init__(self, file_path):
        self.file_path = file_path
        self.hits = BlastHit.generate(self.file_path)
