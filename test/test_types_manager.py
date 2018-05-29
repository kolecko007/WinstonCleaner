#!/usr/bin/env python2.7
import sys, os, unittest
from pathlib import Path

ROOT_PATH = str(Path(os.path.dirname(os.path.realpath(__file__))).parent)
sys.path.append(ROOT_PATH)

from winston.types_manager import TypesManager

class TestTypesManager(unittest.TestCase):
    TYPES_PATH = 'test/test_types.csv'

    def setUp(self):
        with open(os.path.join(ROOT_PATH, self.TYPES_PATH)) as f:
            self.original_content = f.read()

    def tearDown(self):
        with open(os.path.join(ROOT_PATH, self.TYPES_PATH), 'w') as f:
            f.write(self.original_content)

    def test_init(self):
        self.assertIsNone(TypesManager.get_threshold('one', 'two'))
        self.assertIsNone(TypesManager.types_dict)

        TypesManager.load(self.TYPES_PATH)

        result = TypesManager.get_threshold('one', 'two')
        self.assertIsNotNone(result)
        self.assertEqual(result, 0.96)

    def test_order_independency(self):
        TypesManager.load(self.TYPES_PATH)
        print(TypesManager.types_dict)

        thr = 0.96
        self.assertEqual(TypesManager.get_threshold('one', 'two'), thr)
        self.assertEqual(TypesManager.get_threshold('two', 'one'), thr)

        tp = 'CLOSE'
        self.assertEqual(TypesManager.get_type('three', 'four'), tp)
        self.assertEqual(TypesManager.get_type('four', 'three'), tp)


if __name__ == '__main__':
    unittest.main()
