#!/usr/bin/env python2.7
import os, sys, glob
from pathlib import Path
from tqdm import tqdm
from multiprocessing import Pool
from optparse import OptionParser

# TODO: bug in one_vs_one.py line 35, 36
# TODO: fix saving images to the same folder bug

ROOT_PATH = str(Path(os.path.dirname(os.path.realpath(__file__))).parent.parent)
sys.path.append(ROOT_PATH)

from winston.blastab.one_vs_one import OneVsOne
from winston.settings import Settings
from winston.path_resolver import PathResolver

parser = OptionParser()
parser.add_option("-d", "--datasets_path",
                  help="Path for datasets folder with MMETSPXXXX.fas_VS_MMETSPYYYY.fas.blastab files")
parser.add_option("-r", "--results_path", default=ROOT_PATH + '/results/histograms/',
                  help="Path for the histograms folder")
parser.add_option("--left", help="Left organism")
parser.add_option("--right", help="Right organism")

(options, args) = parser.parse_args()


def has_left_right_options():
    return options.left and options.right


def analyze_blastab(path, show=False):
    blastab = OneVsOne(path)
    output_path = '.'.join(path.split('.')[:-2])

    blastab.save_to_file(output_path + '_filtered.blastab')
    blastab.plot(output_path + '.png', show=show)


def main():
    global hist_path
    hist_path = options.results_path
    hist_path = os.path.join(hist_path, '')

    if not os.path.exists(hist_path):
        os.makedirs(hist_path)

    # Formatting input path
    dir_path = options.datasets_path
    dir_path = os.path.join(dir_path, '')

    if has_left_right_options():
        file_path = "%s%s.fas_VS_%s.fas.blastab" % (dir_path, options.left, options.right)
        analyze_blastab(file_path, show=True)
    else:
        files = [fname for fname in glob.glob(dir_path + '*.blastab')]
        files.sort()

        pool = Pool(4)
        for _ in tqdm(pool.imap(analyze_blastab, files), total=len(files)):
            pass


if __name__ == "__main__":
    Settings.load(PathResolver.abs_path_for('config/settings.yml'))
    main()
