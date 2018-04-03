import numpy as np
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt

class HistogramMinimum():
    def __init__(self, histogram, idx, x, y):
        self.histogram = histogram
        self.idx = idx
        self.x = x
        self.y = y

    def prev_bin_value(self):
        index = self.idx-1

        if index <= 0:
            return None

        return self.histogram.bins[self.idx-1]

    def left_border_idx(self):
        for i in range(self.idx, -1, -1):
            if self.histogram.histogram[i] == self.y:
                continue
            else:
                return i

    def left_border_coord(self):
        return self.histogram.x_coord_by_index(self.left_border_idx()+1) - 0.5

    ## rising before the minimum
    def is_correct(self):
        last_id = self.left_border_idx()
        vals = self.histogram.fetch_previous_values(last_id, 2)

        sorted_vals = list(vals)
        sorted_vals.sort()
        sorted_vals.reverse()

        return vals == sorted_vals

    def annotate(self, color='green'):
        plt.annotate("%s\n%s\n%s"%(self.y, self.x, self.idx),
                     xy=(self.x, self.y),
                     xytext=(self.x, self.y),
                     fontsize=6,
                     color='red')
        self.histogram.draw_line(x_coord=self.x, color=color)

class BlastHistogram():
    DEFAULT_HIT_ATTR = 'pident'
    HIST_BINS_COUNT = 50
    DEFAULT_IDENT_THRESHOLD = 95.0
    MINIMUM_HITS_COUNT = 500

    # distanse between MIN and MAX should be equal to bins count
    RANGE_MIN = 51
    RANGE_MAX = 101

    def __init__(self, blastab):
        self.blastab = blastab
        self.hits = self.blastab.hits
        self.values = np.array([getattr(hit, self.DEFAULT_HIT_ATTR) for hit in self.hits])

        rng = (self.RANGE_MIN, self.RANGE_MAX)
        self.histogram, self.bins = np.histogram(self.values, self.HIST_BINS_COUNT, range=rng)

    def len(self):
        return len(self.hits)

    def max_y(self):
        return max(self.histogram)

    def detect_identity_threshold(self):
        if not self.has_enough_hits():
            return self.DEFAULT_IDENT_THRESHOLD

        mins = self.detect_mins()[::-1]
        maxs = self.detect_maxs()[::-1]

        # [m.annotate() for m in mins]
        # mins = [m for m in mins if m.is_correct()]

        dist_threshold = 1
        height_threshold = 20
        correct_mins_count = len([m for m in mins if m.is_correct()])

        current_min = None
        current_min_idx = None
        detected_threshold = None

        for i, m in enumerate(mins):
            if not m.is_correct():
                continue

            if not current_min:
                current_min = m
                current_min_idx = i

                if correct_mins_count == 1:
                    detected_threshold = current_min.left_border_coord()
                    break

                continue

            # 1. next if the height of the max is bad
            if (maxs[i-1].y - mins[current_min_idx].y) > height_threshold:
                detected_threshold = current_min.left_border_coord()
                break

            # 2. next if the distance between mins is bad
            if (current_min.x - m.x) > dist_threshold:
                detected_threshold = current_min.left_border_coord()
                break

            current_min = m
            current_min_idx = i

        if (not detected_threshold) or (detected_threshold and detected_threshold < self.DEFAULT_IDENT_THRESHOLD):
            detected_threshold = self.DEFAULT_IDENT_THRESHOLD

        return detected_threshold

    def plot(self, plt):
        plt.tight_layout()
        plt.hist(self.values, self.HIST_BINS_COUNT, range=(self.RANGE_MIN, self.RANGE_MAX))

        ax = plt.gca()
        ax.set_xticks(np.arange(0, 101, 10))
        ax.set_xticks(np.arange(0, 101, 1), minor=True)

        ax.set_xlabel('pident, %')
        ax.set_ylabel('hits')

        for idx, val in enumerate(self.histogram):
            x = self.x_coord_by_index(idx)
            y = val
            y_loc = self.max_y()+(self.max_y()/7.0) if idx%2 == 0 else self.max_y()+(self.max_y()/6.0)
            plt.annotate(y, xy=(x, y), xytext=(x-1, y_loc), fontsize=4, color='red')

    def detect_mins(self):
        return self.detect_extrema('min')

    def detect_maxs(self):
        return self.detect_extrema('max')

    def detect_extrema(self, target):
        merging_dict = {}
        merged_histogram = []

        # 1. zip the histogram: [1,1,1,2,2,3] -> [1,2,3]
        for idx, val in enumerate(self.histogram):
            if len(merged_histogram) == 0 or val != merged_histogram[-1]:
                merged_histogram.append(val)
                merging_dict[len(merged_histogram)-1] = [idx]

            if idx > 0 and val == self.histogram[idx-1]:
                merging_dict[len(merged_histogram)-1].append(idx)

        # 2. find extrema in a zipped histogram
        target = np.less if target == 'min' else np.greater
        ids = argrelextrema(np.array(merged_histogram), target)[0]

        # 3. unzipping the histogram, calculating coordinates and values
        coords = []
        for idx in ids:
            if len(merging_dict[idx]) == 1:
                x = merging_dict[idx][0] # coord in the source histogram
                m = HistogramMinimum(self, x, self.x_coord_by_index(x), self.histogram[x])
                coords.append(m)
            else:
                x = merging_dict[idx][0]
                avg = np.average([self.x_coord_by_index(i) for i in merging_dict[idx]])
                coords.append(HistogramMinimum(self, x, avg, self.histogram[x]))

        return coords

    def x_coord_by_index(self, idx):
        return np.average([self.bins[idx], self.bins[idx+1]])

    def has_enough_hits(self):
        return self.len() >= self.MINIMUM_HITS_COUNT

    def fetch_previous_values(self, start_idx, count):
        result = []
        current_idx = start_idx

        while len(result) < count:
            if current_idx < 0:
                return result

            current_value = self.histogram[current_idx]

            if len(result) == 0 or result[0] != current_value:
                result.insert(0, current_value)

            current_idx -= 1

        return result


    def draw_line(self, x_coord=None, idx=None, color='red', width=1.0):
        if idx:
            x_coord = self.x_coord_by_index(idx)

        plt.plot((x_coord, x_coord), (0, self.max_y()), 'k-', color=color, linewidth=width)
