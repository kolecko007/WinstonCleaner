import csv

from path_resolver import PathResolver

class TypesManager:
    types_dict = None

    @staticmethod
    def load(file_path=None):
        if not file_path:
            file_path = PathResolver.pair_types_path()

        if TypesManager.types_dict:
            return True

        with open(file_path, 'r') as f:
            result = {}

            for row in csv.reader(f):
                key = TypesManager.detect_key(row[0], row[1])

                result[key] = {
                    'threshold': float(row[2]),
                    'type': row[3]
                    }

        TypesManager.types_dict = result

        return True

    @staticmethod
    def detect_key(left, right):
        return ','.join([left, right])

    @staticmethod
    def get_threshold(left, right):
        return TypesManager.get_type_attr(left, right, 'threshold')

    @staticmethod
    def get_type(left, right):
        return TypesManager.get_type_attr(left, right, 'type')

    @staticmethod
    def get_type_attr(left, right, attr):
        if not TypesManager.types_dict:
            return None
        return TypesManager.types_dict[TypesManager.detect_key(left, right)][attr]

    @staticmethod
    def pair_types_path():
        return PathResolver.output_path_for(PathResolver.TYPES_FILENAME)
