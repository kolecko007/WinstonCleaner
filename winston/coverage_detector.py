import sys, os
import pickle, sqlite3
from tqdm import tqdm
from StringIO import StringIO

from path_resolver import PathResolver

class DatabaseWorker:
    DEFAULT_DB_FILE_NAME = 'database.db'
    TIMEOUT = 10

    IN_MEMORY_DB = None

    @staticmethod
    def load_in_memory_db(db_path=None):
        if not db_path:
            db_path = DatabaseWorker._default_db_path()

        conn = sqlite3.connect(db_path, timeout=DatabaseWorker.TIMEOUT)
        cur = conn.cursor()
        cur.execute("select * from coverage_entries")
        DatabaseWorker.IN_MEMORY_DB = dict(cur.fetchall())

    def __init__(self, db_path=None):
        if not db_path:
            db_path = DatabaseWorker._default_db_path()

        self.db_path = db_path

        if not os.path.exists(db_path):
            open(db_path, 'w').close()

        if not DatabaseWorker.IN_MEMORY_DB:
            self.connection = sqlite3.connect(db_path, timeout=DatabaseWorker.TIMEOUT)
            self.cursor = self.connection.cursor()

    def execute(self, query):
        self.cursor.execute(query)

    def commit(self):
        self.connection.commit()

    @staticmethod
    def _default_db_path():
        return PathResolver.output_path_for(DatabaseWorker.DEFAULT_DB_FILE_NAME)

class CoverageDetector:
    TABLE_NAME = 'coverage_entries'

    def __init__(self, db_path=None):
        self.db_loader = DatabaseWorker(db_path=db_path)

    def create_table(self):
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
        if DatabaseWorker.IN_MEMORY_DB:
            if contig_id in DatabaseWorker.IN_MEMORY_DB:
                return DatabaseWorker.IN_MEMORY_DB[contig_id]
            else:
                return None

        self.db_loader.execute("""
            SELECT * FROM `%s`
            WHERE `contig_id` = '%s';
            """ % (self.TABLE_NAME, contig_id))

        entry = self.db_loader.cursor.fetchone()

        if entry:
            return entry[1]

    def close(self):
        self.db_loader.connection.close()
