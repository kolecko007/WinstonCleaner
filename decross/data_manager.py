import os, glob, itertools, re, subprocess
from pathlib import Path

from settings import Settings
from name_converter import NameConverter
from path_resolver import PathResolver

class DataManager:
    READS_EXT = ['fq', 'fastq']
    READS_FNAME_REGEXP = '^(.+?)_%s\.(' + '|'.join(READS_EXT) + ')$'

    CONTIGS_EXT = ['fa', 'fas', 'fasta']
    CONTIGS_FNAME_REGEXP = '^(.+?)\.(' + '|'.join(CONTIGS_EXT) + ')$'

    def __init__(self):
        self.datasets = []
        self._make_datasets()

    def _make_datasets(self):
        NameConverter.load()

        input_path = str(Path(Settings.decross.paths.input))
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
                if re.match(self.READS_FNAME_REGEXP%1, f):
                    file_group['reads_1'] = f
                elif re.match(self.READS_FNAME_REGEXP%2, f):
                    file_group['reads_2'] = f
                elif re.match(self.CONTIGS_FNAME_REGEXP, f):
                    file_group['contigs'] = f

            if not None in file_group.values():
                dataset = Dataset(file_group['reads_1'],
                                  file_group['reads_2'],
                                  file_group['contigs'])
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
        NameConverter.register(self.external_org_id())

    def prepare(self):
        self._copy_contigs_to_output()
        self._rename_contigs_titles()

    def contigs_input_path(self):
        os.path.join(Settings.decross.paths.input, self.contigs)

    def contigs_output_path(self):
        contigs_file_name = self._get_output_file_name(self.contigs)
        return PathResolver.datasets_output_path(contigs_file_name)

    def reads_input_paths(self):
        return [PathResolver.input_path_for(f) for f in [self.reads_1, self.reads_2]]

    def external_org_id(self):
        return re.search(DataManager.CONTIGS_FNAME_REGEXP, self.contigs).group(1)

    def internal_org_id(self):
        return NameConverter.ext_to_int(self.external_org_id())

    def _copy_contigs_to_output(self):
        old_path = PathResolver.input_path_for(self.contigs)

        new_f_name = self._get_output_file_name(self.contigs)
        new_path = PathResolver.datasets_output_path(new_f_name)

        subprocess.call('cp %s %s' % (old_path, new_path), shell=True)

    def _rename_contigs_titles(self):
        out_file_name = self._get_output_file_name(self.contigs)
        out_file_path = PathResolver.datasets_output_path(out_file_name)

        system_id = self.internal_org_id()
        command = 'sed -i "s#^>\(.*\)\$#>%s_\\1#" %s' % (system_id, out_file_path)
        subprocess.call(command, shell=True)

    def _get_output_file_name(self, input_file_name):
        return input_file_name.replace(self.external_org_id(), self.internal_org_id())

    def _find_file(self, path, external_name, extensions):
        path = os.path.join(path, '')
        paths = [glob.glob('%s%s.%s' % (path, external_name, ext)) for ext in extensions]
        paths = sum(paths, [])
        if len(path) > 0:
            return paths[0]
        else:
            return None


