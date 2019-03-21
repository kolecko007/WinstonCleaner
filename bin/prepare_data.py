#!/usr/bin/env python2.7
import sys, os, logging, glob
from pathlib import Path
from multiprocessing import Pool, Lock
from optparse import OptionParser
from tqdm import tqdm

ROOT_PATH = str(Path(os.path.dirname(os.path.realpath(__file__))).parent)
sys.path.append(ROOT_PATH)

from winston.settings import Settings
from winston.path_resolver import PathResolver
from winston.data_manager import DataManager
from winston.name_converter import NameConverter
from winston.coverage_detector import CoverageDetector
from winston.mapper import Mapper
from winston.blaster import Blaster
from winston.blastab.one_vs_one import OneVsOne

parser = OptionParser(description="Main contamination finder script")
parser.add_option("--config_path", help="Alternative config path")
parser.add_option("--type_detection_only", default=False,
                                           action="store_true",
                                           help="Type detection only")
(options, args) = parser.parse_args()

SETTINGS_PATH = 'settings.yml'

def analyze_blastab(file_path):
    logger.info("Analyzing %s" % file_path)

    blastab = OneVsOne(file_path)
    threshold = blastab.detect_threshold()
    left_name = blastab.left_org_external_name()
    right_name = blastab.right_org_external_name()
    t = blastab.detect_type()
    data = [left_name, right_name, str(threshold), t]

    with outfile_lock:
        with open(PathResolver.pair_types_path(), 'a') as output:
            output.write(",".join(data) + "\n")

    hits_filename = "%s_VS_%s.png" % (left_name, right_name)
    hist_path = PathResolver.output_path_for(PathResolver.HITSOGRAMS_FOLDER, hits_filename)
    blastab.plot(hist_path)

def main():
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)
    global logger
    logger = logging.getLogger()

    Settings.load(options.config_path or SETTINGS_PATH)
    NameConverter.load()

    PathResolver.assure_path_exists(PathResolver.output_path())

    if not options.type_detection_only:
        data_manager = DataManager()

        for i, dataset in enumerate(data_manager.datasets):
            logger.info('Working with dataset %s...' % dataset.external_name)
            dataset.prepare()
            mapper = Mapper(dataset)
            mapper.perform()

            detector = CoverageDetector()
            detector.load_from_bb_tools_rpkm(mapper.pileup_output_path())

        logger.info('BLASTing files')
        blaster = Blaster()
        blaster.perform()

    logger.info('Types detection')
    types_path = PathResolver.output_path_for(PathResolver.TYPES_FILENAME)
    if os.path.exists(types_path):
        os.remove(types_path)

    hist_path = PathResolver.output_path_for(PathResolver.HITSOGRAMS_FOLDER)
    PathResolver.assure_path_exists(hist_path)

    all_vs_all_path = PathResolver.output_path_for(PathResolver.ALL_VS_ALL_FOLDER)
    paths = [fname for fname in glob.glob(all_vs_all_path + '/*.blastab')]
    paths.sort()

    global outfile_lock
    outfile_lock = Lock()

    if Settings.winston.threads.multithreading:
        pool = Pool(int(Settings.winston.threads.count))
        for _ in tqdm(pool.imap(analyze_blastab, paths), total=len(paths)):
            pass
    else:
        for path in paths:
            analyze_blastab(path)

if __name__ == '__main__':
    main()
