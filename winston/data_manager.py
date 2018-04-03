import os, glob, itertools, re, subprocess
from pathlib import Path

from settings import Settings
from name_converter import NameConverter
from path_resolver import PathResolver
from Bio import SeqIO

class DataManager:
    READS_EXT = ['fq', 'fastq']
    READS_FNAME_REGEXP = '^(.+?)_%s\.(' + '|'.join(READS_EXT) + ')$'
    LEFT_READS_FNAME_REGEXP = re.compile(READS_FNAME_REGEXP % 1)
    RIGHT_READS_FNAME_REGEXP = re.compile(READS_FNAME_REGEXP % 2)

    CONTIGS_EXT = ['fa', 'fas', 'fasta']
    CONTIGS_FNAME_REGEXP = re.compile('^(.+?)\.(' + '|'.join(CONTIGS_EXT) + ')$')

    def __init__(self):
        self.datasets = []
        self._make_datasets()

    def _make_datasets(self):
        NameConverter.load()

        input_path = PathResolver.input_path()
        files = [n.split('/')[-1] for n in glob.glob(input_path + '/*')]
        files.sort()

        groups = {}
        for gid, els in itertools.groupby(files, self._extract_external_name):
            groups[gid] = list(els)

        for group, files in groups.iteritems():
            if not group or len(files) < 3:
                continue

            file_group = { 'reads_1': None, 'reads_2': None, 'contigs': None }

            for f in files:
                if re.match(DataManager.LEFT_READS_FNAME_REGEXP, f):
                    file_group['reads_1'] = f
                elif re.match(DataManager.RIGHT_READS_FNAME_REGEXP, f):
                    file_group['reads_2'] = f
                elif re.match(DataManager.CONTIGS_FNAME_REGEXP, f):
                    file_group['contigs'] = f


            if not None in file_group.values():
                ext_name = self._extract_external_name(file_group['contigs'])
                dataset = Dataset(ext_name)

                self.datasets.append(dataset)

        return self.datasets


    def _extract_external_name(self, file_name):
        result = re.search('(.+?)(_[12])?\..+', file_name)
        if result:
            return result.group(1)
        else:
            return None

class Dataset:
    """ Logic of input files for the organism
    """

    def __init__(self, external_name):
        if not os.path.isdir(PathResolver.datasets_output_path()):
            os.makedirs(PathResolver.datasets_output_path())

        self.external_name = external_name
        self.internal_name = self.detect_internal_name()
        NameConverter.register(self.external_name)

        self._contigs_count = None

    def prepare(self):
        self._copy_contigs_to_output()
        self._rename_contigs_titles()

    # Input paths

    def contigs_input_path(self):
        ext = DataManager.CONTIGS_EXT
        return self._find_dataset_file(self.external_name, ext)

    def reads_input_paths(self):
        file_names = ["%s_%s" % (self.external_name, i) for i in [1, 2]]
        ext = DataManager.READS_EXT
        return [self._find_dataset_file(fn, ext) for fn in file_names]

    # Output paths

    def contigs_output_path(self):
        file_name = os.path.basename(self.contigs_input_path())
        return PathResolver.datasets_output_path(file_name)

    def reads_output_paths(self):
        input_names = [os.path.basename(p) for p in self.reads_input_paths()]
        return [PathResolver.datasets_output_path(p) for p in input_names]

    def detect_internal_name(self):
        return NameConverter.ext_to_int(self.external_name)

    def contigs_count(self):
        if not self._contigs_count:
            path = self.contigs_input_path()
            cnt = subprocess.check_output('grep -c ">" %s' % path, shell=True)
            self._contigs_count = int(cnt)

        return self._contigs_count

    def _copy_contigs_to_output(self):
        old_path = self.contigs_input_path()
        new_path = self.contigs_output_path()

        subprocess.call('cp %s %s' % (old_path, new_path), shell=True)

    def _rename_contigs_titles(self):
        out_file_path = self.contigs_output_path()
        tmp_file_path = out_file_path + '.tmp'
        system_id = self.internal_name

        with open(out_file_path) as out_f, open(tmp_file_path, 'w') as tmp_f:
            for record in SeqIO.parse(out_f, 'fasta'):
                new_id = record.id
                new_id = re.split('\s+', new_id)[0]
                new_id = "%s_%s" % (system_id, new_id)
                record.id = new_id
                record.description = ''
                SeqIO.write(record, tmp_f, 'fasta')

        subprocess.call('mv %s %s' % (tmp_file_path, out_file_path), shell=True)

    def _get_output_path(self, input_path):
        file_name = os.path.basename(input_path)
        file_name = input_file_name.replace(file_name, self.internal_name)
        return PathResolver.datasets_output_path(file_name)

    def _find_dataset_file(self, file_name, extensions):
        path = PathResolver.input_path()
        paths = [glob.glob('%s/%s.%s' % (path, file_name, ext)) for ext in extensions]
        paths = sum(paths, [])
        if len(paths) > 0:
            return paths[0]
        else:
            return None


