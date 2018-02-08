import sys, os
import pickle, sqlite3
from tqdm import tqdm
from StringIO import StringIO

from path_resolver import PathResolver

class DatabaseWorker:
    DEFAULT_DB_FILE_NAME = 'database.db'
    TIMEOUT = 10

    IN_MEMORY_CONNECTION = None

    def __init__(self, db_path=None, in_memory=False):
        if not db_path:
            db_path = self._default_db_path()

        self.db_path = db_path

        if not os.path.exists(db_path):
            open(db_path, 'w').close()

        if in_memory:
            if not DatabaseWorker.IN_MEMORY_CONNECTION:
                print 'init connection!!!!'
                print
                con = sqlite3.connect(db_path, timeout=self.TIMEOUT)
                tempfile = StringIO()
                for line in con.iterdump():
                    tempfile.write('%s\n' % line)
                con.close()
                tempfile.seek(0)

                DatabaseWorker.IN_MEMORY_CONNECTION = sqlite3.connect(":memory:")
                DatabaseWorker.IN_MEMORY_CONNECTION.cursor().executescript(tempfile.read())
                DatabaseWorker.IN_MEMORY_CONNECTION.commit()
                DatabaseWorker.IN_MEMORY_CONNECTION.row_factory = sqlite3.Row

            self.connection = DatabaseWorker.IN_MEMORY_CONNECTION
            self.cursor = DatabaseWorker.IN_MEMORY_CONNECTION.cursor()
        else:
            self.connection = sqlite3.connect(db_path, timeout=self.TIMEOUT)
            self.cursor = self.connection.cursor()

    def execute(self, query):
        self.cursor.execute(query)

    def commit(self):
        self.connection.commit()

    def _default_db_path(self):
        return PathResolver.output_path_for(self.DEFAULT_DB_FILE_NAME)

class CoverageDetector:
    TABLE_NAME = 'coverage_entries'

    def __init__(self, db_path=None, in_memory=False):
        self.db_loader = DatabaseWorker(db_path=db_path, in_memory=in_memory)

    def create_table(self):
        db_path = self.db_loader.db_path

        if os.path.exists(db_path):
            os.remove(db_path)
            self.db_loader = DatabaseWorker(db_path=db_path)

        self.db_loader.execute("""
            create table if not exists `%s` (
            `contig_id` char[256] NOT NULL,
            `kmer` float NOT NULL,
            PRIMARY KEY (`contig_id`))
            """ % self.TABLE_NAME)

    def insert(self, contig_id, kmer):
        self.db_loader.execute("""
            INSERT INTO `%s` (`contig_id`, `kmer`)
            VALUES ('%s', '%s');
            """ % (self.TABLE_NAME, contig_id, kmer))

    def load_from_pickle(self, file_path):
        items = pickle.load(open(file_path, 'r'))
        self.load_from_hash(items)

    def load_from_csv(self, file_path):
        items = {}
        file = open(file_path, 'r')
        for line in file.readlines():
            splitted = line.split(',')
            items[splitted[0]] = [float(splitted[1])]
        self.load_from_hash(items)

    def load_from_bb_tools_rpkm(self, file_path):
        items = {}

        with open(file_path, 'r') as f:
            for line in f.readlines():
                if line.startswith('#'):
                    continue
                splitted = line.split('\t')
                items[splitted[0].split(' ')[0]] = [float(splitted[5])]

        self.load_from_hash(items)

    def load_from_hash(self, items):
        self.create_table()
        commit_items_count = 5000
        cnt = 0

        for name, value in tqdm(items.iteritems(), total=len(items)):
            self.insert(name, value[0])

            if cnt >= commit_items_count:
                cnt = 0
                self.db_loader.commit()

            cnt += 1

        self.db_loader.commit()

    def kmer_by_contig_id(self, contig_id):
        self.db_loader.execute("""
            SELECT * FROM `%s`
            WHERE `contig_id` = '%s';
            """ % (self.TABLE_NAME, contig_id))

        entry = self.db_loader.cursor.fetchone()

        if entry:
            return entry[1]

    def close(self):
        self.db_loader.connection.close()
