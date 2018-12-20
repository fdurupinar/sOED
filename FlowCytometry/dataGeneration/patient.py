import numpy as np
from util.staticMethods import StaticMethods
import matplotlib.pyplot as plt


class Patient:

    def __init__(self, is_cancer, marker_list, marker_cnt, cell_cnt):

        self.is_cancer = is_cancer
        self.marker_list = marker_list
        self.marker_cnt = {}
        self.cell_cnt = cell_cnt

        for i in range(0, len(marker_cnt)):
            key = StaticMethods.get_ab_key([marker_list[i]])
            self.marker_cnt[key] = marker_cnt[i]



    def get_marker_count(self, ab_list):
        """
        For all the cells, count the ones with all the markers in ab_seq and return their ratios
        :param ab_list: Is in format [A,B,C] : means antibodies A, B, and C are present
        :param to_be_recorded:
        :return:
        """

        key = StaticMethods.get_ab_key(ab_list)

        return self.marker_cnt[key]
