import os, subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from ..path_resolver import PathResolver
from ..name_converter import NameConverter
from ..blast_hit import BlastHit
from ..blast_histogram import BlastHistogram
from ..data_manager import Dataset
from blastab import Blastab

class OneVsOne(Blastab):
    PLOT_DPI = 200
    PP_THRESHOLD = 40 # %
    PP_CONT_THRESHOLD = 10 # %
    IDENTICAL_SPECIES_PERCENTAGE = 5 # %
    IDENTICAL_HITS_THRESHOLD = 5 # %, to be identical, hits count must be >= IDENTICAL_HITS_THRESHOLD of left orgranism dataset count

    TYPES = { 1: 'NO', 2: 'REGULAR', 3: 'CLOSE', 4: 'LEFT_EATS_RIGHT', 5: 'RIGHT_EATS_LEFT' }

    def __init__(self, file_path):
        super(OneVsOne, self).__init__(file_path)

        self.merge_hits()
        self.filter_hits()

        self.histogram = BlastHistogram(self)
        self.dataset_dir = PathResolver.datasets_output_path()

        self.threshold = None
        self.type = None

        self.left_dataset = Dataset(self.left_org_external_name())
        self.right_dataset = Dataset(self.right_org_external_name())

    def detect_threshold(self):
        if self.threshold:
            return self.threshold

        if self.is_close():
            return BlastHistogram.DEFAULT_IDENT_THRESHOLD

        self.threshold = self.histogram.detect_identity_threshold()

        return self.threshold

    def detect_type(self):
        if self.type:
            return self.type

        if self.is_possibly_identical():
            self.type = self.TYPES[3] # CLOSE
        else:
            self.type = self.TYPES[2] # REGULAR

        return self.type

    def plot(self, output_path, show=False):
        left_name = self.left_org_external_name()
        right_name = self.right_org_external_name()

        title = left_name + " VS " + right_name + " (threshold: " + str(self.threshold) + ", " + str(len(self.histogram.hits)) + " hits)"
        plt.title(title, color='red', fontsize=8)

        self.histogram.plot(plt)

        if self.threshold:
            self.histogram.draw_line(x_coord=self.threshold)

        plt.savefig(output_path, dpi=self.PLOT_DPI)

        if show:
            plt.show()

        plt.close()

    def save_to_file(self, path):
        with open(path, 'w') as f:
            for hit in self.hits:
                f.write(hit.raw)

    def left_org_external_name(self):
        return self.org_external_names()[0]

    def right_org_external_name(self):
        return self.org_external_names()[1]

    def org_external_names(self):
        return os.path.splitext(os.path.basename(self.file_path))[0].split('_VS_')

    def left_dataset_count(self):
        return self.left_dataset.contigs_count()

    def right_dataset_count(self):
        return self.right_dataset.contigs_count()

    def merge_hits(self):
        ''' Keeps only best hits '''
        name_hash = {}

        for hit in self.hits:
            if hit.sseqid not in name_hash:
                name_hash[hit.sseqid] = []
            name_hash[hit.sseqid].append(hit)

        merged = []

        for name, array in name_hash.iteritems():
            if len(array) == 0:
                raise Exception('Error: array should not be empty')

            if len(array) == 1:
                merged.append(array[0])
                continue

            if len(array) > 1:
                best = None
                for hit in array:
                    if not best:
                        best = hit
                    else:
                        if hit.length > best.length:
                            best = hit

                merged.append(best)

        self.hits = merged
        return self.hits

    def filter_hits(self):
        self.hits = [x for x in self.hits if x.is_reliable()]
        return self.hits

    def is_regular(self):
        return self.detect_type() == self.TYPES[2]

    def is_close(self):
        return self.detect_type() == self.TYPES[3]

    def is_absolutely_identical(self):
        return self.left_org_external_name() == self.right_org_external_name()

    def is_possibly_identical(self):
        left = sum(self.histogram.histogram[:-5])
        right = sum(self.histogram.histogram[-5:])

        total_count = float(left+right)
        if total_count > 0:
            percentage = float(left)/total_count
        else:
            percentage = 1.0

        identical_by_histogram = percentage <= self.IDENTICAL_SPECIES_PERCENTAGE/100.0

        contigs_count = self.left_dataset_count()
        enough_hits = len(self.hits) >= contigs_count*(self.IDENTICAL_HITS_THRESHOLD/100.0)

        return identical_by_histogram and enough_hits


