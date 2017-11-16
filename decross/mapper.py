import os, subprocess

from settings import Settings
from path_resolver import PathResolver

class Mapper:
    MAPPINGS_FOLDER = 'coverage'

    def __init__(self, dataset):
        self.dataset = dataset

        folder_path = PathResolver.output_path_for(self.MAPPINGS_FOLDER)
        if not os.path.isdir(folder_path):
            os.makedirs(folder_path)

    def perform(self):
        self._map()
        self._pileup()

    def sam_file_path(self):
        file_name = '%s.sam' % self.dataset.external_name
        return self._mappings_folder_path(file_name)

    def pileup_output_path(self):
        file_name = '%s_pileup.txt' % self.dataset.external_name
        return self._mappings_folder_path(file_name)

    def _map(self):
        try:
            exe = Settings.decross.paths.tools.bbmap_sh
        except AttributeError:
            exe = 'bbmap.sh'

        command = exe + ' in=%s in2=%s ref=%s nodisk out=%s'
        in_1, in_2 = self.dataset.reads_input_paths()
        ref = self.dataset.contigs_output_path()
        out = self.sam_file_path()
        subprocess.call(command % (in_1, in_2, ref, out), shell=True)

    def _pileup(self):
        try:
            exe = Settings.decross.paths.tools.pileup_sh
        except AttributeError:
            exe = 'pileup.sh'

        command = exe + " in=%s ref=%s rpkm=%s overwrite=true"
        sam_path = self.sam_file_path()
        ref = self.dataset.contigs_output_path()
        rpkm = self.pileup_output_path()
        subprocess.call(command % (sam_path, ref, rpkm), shell=True)

    def _mappings_folder_path(self, file_name):
        return PathResolver.output_path_for(self.MAPPINGS_FOLDER, file_name)
