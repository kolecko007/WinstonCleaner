#!/usr/bin/python2.7
import sys, os, progressbar, logging, glob
from pathlib import Path
from multiprocessing import Pool, Lock
from optparse import OptionParser
from tqdm import tqdm

ROOT_PATH = str(Path(os.path.dirname(os.path.realpath(__file__))).parent)
sys.path.append(ROOT_PATH)

from decross.settings import Settings
from decross.path_resolver import PathResolver
from decross.data_manager import DataManager
from decross.name_converter import NameConverter
from decross.coverage_detector import CoverageDetector
from decross.mapper import Mapper
from decross.blaster import Blaster
from decross.blastab.one_vs_one import OneVsOne

parser = OptionParser(description="Main contamination finder script")
parser.add_option("--config_path", help="Alternative config path")
parser.add_option("--type_detection_only", default=False,
                                           action="store_true",
                                           help="Type detection only")
(options, args) = parser.parse_args()

SETTINGS_PATH = 'config/settings.yml'

def analyze_blastab(file_path):
    logger.info("Analyzing %s" % file_path)

    blastab = OneVsOne(file_path)
    threshold = blastab.detect_threshold()
    left_name = blastab.left_org_outer_name()
    right_name = blastab.right_org_outer_name()
    t = blastab.detect_type()
    data = [left_name, right_name, str(threshold), t]

    with oufile_lock:
        with open(PathResolver.pair_types_path(), 'a') as output:
          output.write(",".join(data) + "\n")

    hits_filename = "%s_vs_%s.png" % (left_name, right_name)
    hist_path = PathResolver.output_path_for(PathResolver.HITSOGRAMS_FOLDER, hits_filename)
    blastab.plot(hist_path)

def main():
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)
    global logger
    logger = logging.getLogger()

    Settings.load(options.config_path or PathResolver.abs_path_for(SETTINGS_PATH))
    NameConverter.load()

    output_path = Settings.decross.paths.output
    PathResolver.assure_path_exists(output_path)

    if not options.type_detection_only:
        data_manager = DataManager()

        for i, dataset in enumerate(data_manager.datasets):
            logger.info('Working with dataset %s...' % dataset.external_org_id())
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

    global oufile_lock
    oufile_lock = Lock()

    if Settings.decross.threads.multithreading:
        pool = Pool(int(Settings.decross.threads.count))
        for _ in tqdm(pool.imap(analyze_blastab, paths), total=len(paths)):
            pass
    else:
        for path in paths:
            analyze_blastab(path)

if __name__ == '__main__':
    main()
