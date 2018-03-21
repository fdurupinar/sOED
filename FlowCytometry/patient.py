import numpy as np
from staticMethods import StaticMethods
import matplotlib.pyplot as plt


class Patient:

    def __init__(self, cell_cnt, ab_cnt, markers_list, mu_list, sigma_list):

        self.cells = np.zeros([cell_cnt, ab_cnt])  # each cell has ab_cnt bins, which can be filled with markers or not
        self.ab_cnt = ab_cnt
        # self.ab_percentage_arr = np.full((ab_cnt), -1, dtype=np.int)
        self.is_marker_arr = np.full((ab_cnt), False, dtype=np.bool)
        self.mu_list = mu_list
        self.sigma_list = sigma_list
        self.markers_list = markers_list

        # self._compute_ab_percentages()
        # self._assign_all_antibodies()

        self._fill_in_cells()
        self.marker_ratios = {}

    def _fill_in_cells(self):
        """
        Fill in the percentages array for each antibody
        :return:
        """

        # first fill in the percentages for markers
        for i in range(len(self.markers_list)):
            prob = np.random.normal(self.mu_list[i], self.sigma_list[i])
            for m in self.markers_list[i]:
                # assign cells of that antibody group to 1 with that percentage
                for j in range(len(self.cells)):
                    val = np.random.uniform(0, 1)
                    if val <= prob:
                        self.cells[j][m] = 1
                self.is_marker_arr[m] = True

        # then fill in the rest of the antibodies drawing their ratios from a uniform distribution
        for i in range(self.ab_cnt):  # for each protein
            if not self.is_marker_arr[i]:  # not in the marker group
                # TODO uniform dist
                # prob = np.random.uniform(0, 1)
                # prob = 0.5
                prob = 0
                for j in range(len(self.cells)):
                    val = np.random.uniform(0, 1)
                    if val <= prob:
                        self.cells[j][i] = 1

    def get_marker_ratio(self, ab_list, to_be_recorded):
        """
        For all the cells, count the ones with all the markers in ab_seq and return their ratios
        :param ab_list: Is in format [A,B,C] : means antibodies A, B, and C are present
        :param to_be_recorded:
        :return:
        """

        key = StaticMethods.get_ab_key(ab_list)
        if key in self.marker_ratios:
            return self.marker_ratios[key]

        containing_cell_cnt = 0
        for c in self.cells:
            cnt_pres = 0
            for ab in ab_list:
                if int(c[ab]) == 1:
                    cnt_pres += 1
            if cnt_pres >= len(ab_list):
                    containing_cell_cnt += 1

        ratio = float(containing_cell_cnt) / len(self.cells)

        if to_be_recorded:
            self.marker_ratios[key] = ratio
        return ratio

    def plot_patient(self):

        x_data = []
        y_data = []

        for i in range(self.ab_cnt):

            x_data.append(i)
            y_data.append(self.get_marker_ratio([i], False))

        plt.plot(x_data, y_data)

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
