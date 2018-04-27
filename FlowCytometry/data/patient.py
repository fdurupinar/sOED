import numpy as np
from util.staticMethods import StaticMethods
import matplotlib.pyplot as plt


class Patient:

    def __init__(self, is_cancer, cell_cnt, ab_cnt, markers_list, mu_list, sigma_list):

        self.cells = np.zeros([cell_cnt, ab_cnt])  # each cell has ab_cnt bins, which can be filled with markers or not
        self.cell_cnt = cell_cnt
        self.ab_cnt = ab_cnt
        self.is_marker_arr = np.full((ab_cnt), False, dtype=np.bool)
        self.mu_list = mu_list
        self.sigma_list = sigma_list
        self.markers_list = markers_list
        self.is_cancer = is_cancer

        self._fill_in_cells()
        self.marker_cnt = {}

        # self.plot_patient()

    def _fill_in_cells(self):
        """
        Fill in the percentages array for antibody groups
        :return:
        """

        # first in the  antibodies drawing their ratios from a uniform distribution
        for i in range(self.ab_cnt):  # for each protein
                # TODO uniform dist
            # prob = np.random.uniform(0, 1)
                # prob = 0.5
            prob = 0
            for j in range(self.cell_cnt):
                val = np.random.uniform(0, 1)
                if val <= prob:
                    self.cells[j][i] = 1

        # then create a bias towards the markers

        # make all the markers true

        for i in range(len(self.markers_list)):
            for m in self.markers_list[i]:
                # assign cells of that antibody group to val with that percentage
                for j in range(self.cell_cnt):
                    prob = np.random.normal(self.mu_list[i], self.sigma_list[i])
                    val = np.random.uniform(0, 1)
                    if val <= prob:
                        self.cells[j][m] = 1
                self.is_marker_arr[m] = True

        # remove that bias in nc patients by picking a random marker and removing it
        if not self.is_cancer:
            for i in range(len(self.markers_list)):
                # assign cells of that antibody group to val with that percentage
                for j in range(self.cell_cnt):
                    # make one of the markers in the list false for all the cells
                    ind = np.random.randint(len(self.markers_list[i]))
                    m = self.markers_list[i][ind]

                    self.cells[j][m] = 0

        # print self.cells

    def get_marker_count(self, ab_list, to_be_recorded):
        """
        For all the cells, count the ones with all the markers in ab_seq and return their ratios
        :param ab_list: Is in format [A,B,C] : means antibodies A, B, and C are present
        :param to_be_recorded:
        :return:
        """

        key = StaticMethods.get_ab_key(ab_list)
        if key in self.marker_cnt:
            return self.marker_cnt[key]

        containing_cell_cnt = 0
        for c in self.cells:
            cnt_pres = 0
            for ab in ab_list:
                if int(c[ab]) == 1:
                    cnt_pres += 1
            if cnt_pres >= len(ab_list):
                    containing_cell_cnt += 1

        cnt = float(containing_cell_cnt)

        if to_be_recorded:
            self.marker_cnt[key] = cnt
        return cnt

    def plot_patient(self):

        x_data = np.arange(self.ab_cnt)
        y_data = []
        for i in range(self.ab_cnt):
            y_data.append(self.get_marker_count([i], False))

        plt.bar(x_data, y_data)

        plt.xlabel('Antibody index')
        plt.ylabel('Marker ratio')
        plt.grid(True)

        plt.show()


# pt = Patient(100, 240, [[1, 2, 3, 4]], [0.8], [0.1])
# pt.plot_patient()
# print pt.cells
# print pt.get_marker_ratio([1,2], [])

# for c in pt.cells:
#     print c
