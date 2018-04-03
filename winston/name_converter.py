import os, random

from settings import Settings
from path_resolver import PathResolver

class NameConverter:
    DICT_FILENAME = 'system_names.csv'
    ext_to_int_data = {}
    int_to_ext_data = {}

    @staticmethod
    def ext_to_int(ext_name):
        if ext_name not in NameConverter.ext_to_int_data:
            NameConverter.register(ext_name)

        return NameConverter.ext_to_int_data[ext_name]

    @staticmethod
    def int_to_ext(int_name):
        if int_name not in NameConverter.int_to_ext_data:
            return None

        return NameConverter.int_to_ext_data[int_name]

    @staticmethod
    def register(ext_name):
        NameConverter.load()

        if ext_name in NameConverter.ext_to_int_data:
            return False

        int_name = "%032x" % random.getrandbits(128)

        with open(NameConverter._dict_file_path(), 'a') as f:
            f.write("%s,%s\n" % (ext_name, int_name))
            NameConverter.ext_to_int_data[ext_name] = int_name
            NameConverter.int_to_ext_data[int_name] = ext_name

        return True

    @staticmethod
    def load():
        NameConverter.int_to_ext_data = {}
        NameConverter.ext_to_int_data = {}

        if not os.path.exists(NameConverter._dict_file_path()):
            return False

        with open(NameConverter._dict_file_path(), 'r') as f:
            for line in f.readlines():
                vals = line.strip().split(',')
                NameConverter.ext_to_int_data[vals[0]] = vals[1]
                NameConverter.int_to_ext_data[vals[1]] = vals[0]

        return True

    @staticmethod
    def _dict_file_path():
        return PathResolver.output_path_for(NameConverter.DICT_FILENAME)
