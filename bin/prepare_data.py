#!/usr/bin/python2.7

import sys, os, progressbar, logging
from pathlib import Path

ROOT_PATH = str(Path(os.path.dirname(os.path.realpath(__file__))).parent)
sys.path.append(ROOT_PATH)

from decross.settings import Settings
from decross.path_resolver import PathResolver
from decross.input_manager import InputManager
from decross.coverage_detector import CoverageDetector
from decross.mapper import Mapper

SETTINGS_PATH = 'config/settings.yml'

def process_dataset(dataset):
    dataset.prepare()

    mapper = Mapper(dataset)
    mapper.perform()

    detector = CoverageDetector()
    detector.load_from_bb_tools_rpkm(mapper.pileup_output_path())

def main():
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)
    logger = logging.getLogger()

    Settings.load(PathResolver.abs_path_for(SETTINGS_PATH))

    output_path = Settings.decross.paths.output
    PathResolver.assure_path_exists(output_path)

    input_manager = InputManager()

    for i, dataset in enumerate(input_manager.datasets):
        logger.info('Working with dataset %s...' % dataset.external_org_id())
        process_dataset(dataset)

    logger.info('BLASTing files')
    blaster = Blaster()
    blaster.perform()

if __name__ == '__main__':
    main()
