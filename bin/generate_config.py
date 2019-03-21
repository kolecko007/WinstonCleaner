#!/usr/bin/env python2.7
import os
import sys
import shutil
from pathlib import Path
from optparse import OptionParser

ROOT_PATH = str(Path(os.path.dirname(os.path.realpath(__file__))).parent)
sys.path.append(ROOT_PATH)

from winston.path_resolver import PathResolver

SETTINGS_PATH = 'winston/config/settings.yml.default'


def parse_args():
    parser = OptionParser(description="Generate default config")
    parser.add_option("--output_folder", help="Where to put new settings.yml (current folder by default)")
    (options, args) = parser.parse_args()
    return options


def main():
    options = parse_args()
    dummy_path = PathResolver.abs_path_for(SETTINGS_PATH)
    output_path = "settings.yml"

    if options.output_folder:
        output_path = os.path.join(options.output_folder, output_path)

    shutil.copyfile(dummy_path, output_path)

    print("New dummy config has succesfully copied to %s" % output_path)


if __name__ == '__main__':
    main()
