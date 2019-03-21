import os, subprocess

from settings import Settings
from path_resolver import PathResolver

class Mapper:
    MAPPINGS_FOLDER = 'coverage'
    DB_FOLDER = 'db'

    def __init__(self, dataset):
        self.dataset = dataset

        folder_path = PathResolver.output_path_for(self.MAPPINGS_FOLDER)
        PathResolver.assure_path_exists(folder_path)

        db_path = PathResolver.output_path_for(self.MAPPINGS_FOLDER, self.DB_FOLDER)
        PathResolver.assure_path_exists(db_path)

    def perform(self):
        self._map()
        self._pileup()

    def sam_file_path(self):
        file_name = '%s.sam' % self.dataset.external_name
        paths = [self.MAPPINGS_FOLDER, file_name]
        return PathResolver.output_path_for(*paths)

    def pileup_output_path(self):
        file_name = '%s_pileup.txt' % self.dataset.external_name
        return PathResolver.output_path_for(self.MAPPINGS_FOLDER, file_name)

    def db_path(self):
        return PathResolver.output_path_for(self.MAPPINGS_FOLDER, self.DB_FOLDER)

    def _map(self):
        try:
            bt_exe = Settings.winston.paths.tools.bowtie2
        except AttributeError:
            bt_exe = 'bowtie2'

        try:
            build_exe = Settings.winston.paths.tools.bowtie2_build
        except AttributeError:
            build_exe = 'bowtie2-build'

        in_1, in_2 = self.dataset.reads_input_paths()
        ref_path = self.dataset.contigs_output_path()
        db_path = self.db_path()
        out_path = self.sam_file_path()

        command = build_exe + ' %s %s'
        command = command % (ref_path, db_path)
        subprocess.call(command, shell=True)

        threads_cnt = Settings.winston.tools.bowtie2.threads
        command = bt_exe + ' --very-sensitive -p %s -x %s -1 %s -2 %s -S %s'
        command = command % (threads_cnt, db_path, in_1, in_2, out_path)

        subprocess.call(command, shell=True)

    def _pileup(self):
        try:
            exe = Settings.winston.paths.tools.pileup_sh
        except AttributeError:
            exe = 'pileup.sh'

        command = exe + " in=%s ref=%s rpkm=%s overwrite=true"
        sam_path = self.sam_file_path()
        ref = self.dataset.contigs_output_path()
        rpkm = self.pileup_output_path()
        subprocess.call(command % (sam_path, ref, rpkm), shell=True)
