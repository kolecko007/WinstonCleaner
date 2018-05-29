#!/usr/bin/env python2.7
import unittest

ROOT_PATH = str(Path(os.path.dirname(os.path.realpath(__file__))).parent)
sys.path.append(ROOT_PATH)

class DataManagerTest(unittest.TestCase):
    def setUp(self):
        print 'up'
        pass

    def tearDown(self):
        print 'down'
        pass

    def test_foo(self):
        print 'foo'
        pass

    def test_bar(self):
        print 'bar'
        pass

if __name__ == '__main__':
    unittest.main()
