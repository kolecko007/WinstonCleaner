import os, glob, re, itertools, subprocess

from path_resolver import PathResolver
from data_manager import DataManager
from settings import Settings

class Blaster:
    COLUMNS = 'qseqid sseqid pident length qlen gaps slen qstart qend qcovs qcovhsp'
    BLAST_RESULT_EXT = 'blastab'

    def __init__(self):
        pass

    def perform(self):
        PathResolver.assure_path_exists(self._db_folder_path())

        all_vs_all_path = PathResolver.output_path_for(PathResolver.ALL_VS_ALL_FOLDER)
        PathResolver.assure_path_exists(all_vs_all_path)

        one_vs_all_path = PathResolver.output_path_for(PathResolver.ONE_VS_ALL_FOLDER)
        PathResolver.assure_path_exists(one_vs_all_path)

        self._make_databases()
        self._perform_blasts()
        self._make_one_vs_all_files()

    def get_all_vs_all_pairs(self):
        product = list(itertools.product(self._contigs_files(), repeat=2))
        return [x for x in product if x[0] != x[1]]

    def _contigs_files(self):
        files = glob.glob(PathResolver.datasets_output_path() + '/*')
        return [e for e in files if re.match(DataManager.CONTIGS_FNAME_REGEXP, e)]

    def _make_databases(self):
        for f_path in self._contigs_files():
            self._make_database(f_path)

    def _make_database(self, file_path):
        command = 'makeblastdb -in %s -parse_seqids -dbtype nucl -out %s'

        org_name = os.path.splitext(file_path)[0].split('/')[-1]
        db_path = self._db_folder_path(org_name)

        subprocess.call(command % (file_path, db_path), shell=True)

    def _perform_blasts(self):
        for left, right in self.get_all_vs_all_pairs():
            self._perform_blast(left, right)

    def _perform_blast(self, left_path, right_path):
        command = 'blastn -query %s -db %s -out %s -outfmt "6 %s" -num_threads %s'

        right_org_name = os.path.splitext(os.path.basename(right_path))[0]
        left_org_name = os.path.splitext(os.path.basename(left_path))[0]
        db_path = self._db_folder_path(right_org_name)

        outfile_name = '%s_VS_%s.%s' % (left_org_name, right_org_name, self.BLAST_RESULT_EXT)
        output_path = PathResolver.output_path_for(PathResolver.ALL_VS_ALL_FOLDER, outfile_name)

        threads_cnt = Settings.winston.tools.blast.threads
        command = command % (left_path, db_path, output_path, self.COLUMNS, threads_cnt)

        subprocess.call(command, shell=True)

    def _make_one_vs_all_files(self):
        all_vs_all_path = PathResolver.output_path_for(PathResolver.ALL_VS_ALL_FOLDER)

        files = {}
        for f_path in glob.glob(all_vs_all_path + '/*.%s' % self.BLAST_RESULT_EXT):
            left_org_id = os.path.splitext(os.path.basename(f_path))[0].split('_VS_')[0]

            if left_org_id not in files:
                files[left_org_id] = []
            files[left_org_id].append(f_path)

        for org_id, files in files.iteritems():
            org_file_name = org_id + '.' + self.BLAST_RESULT_EXT
            path = PathResolver.output_path_for(PathResolver.ONE_VS_ALL_FOLDER, org_file_name)

            if os.path.exists(path):
                os.remove(path)

            for f_path in files:
                subprocess.call('cat %s >> %s' % (f_path, path), shell=True)

    def _db_folder_path(self, *inner_path):
        return PathResolver.output_path_for(PathResolver.BLAST_DB_FOLDER, *inner_path)
