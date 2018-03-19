import numpy as np
from staticMethods import StaticMethods
import matplotlib.pyplot as plt


class Patient:

    def __init__(self, cell_cnt, ab_cnt, markers, mu, sigma):

        self.cells = np.zeros([cell_cnt, ab_cnt])  # each cell has ab_cnt bins, which can be filled with markers or not
        self.ab_cnt = ab_cnt
        self.ab_percentage_arr = np.zeros([ab_cnt])
        self.mu = mu
        self.sigma = sigma
        self.markers = markers

        self._compute_ab_percentages()
        self._assign_all_antibodies()

        self.marker_ratios = {}

    def _compute_ab_percentages(self):
        """
        Fill in the percentages array for each antibody
        :return:
        """

        for i in range(self.ab_cnt):  # for each protein
            if i in self.markers:
                self.ab_percentage_arr[i] = np.random.normal(self.mu,
                                                             self.sigma)  # draw marker probability from a normal distribution
            else:  # draw the rest of the probabilities from a uniform distribution
                #TODO try uniform dist
                # self.ab_percentage_arr[i] = np.random.uniform(0, 1)
                self.ab_percentage_arr[i] = 0

            if self.ab_percentage_arr[i] < 0:
                self.ab_percentage_arr[i] = 0
            if self.ab_percentage_arr[i] > 1:
                self.ab_percentage_arr[i] = 1

    def _assign_ab(self, ab, prob):
        """
        Assign an antibody with index ab with a probability prob
        :param ab: Antibody id 
        :param prob: Probability between 0, 1 
        :return: 
        """

        for i in range(len(self.cells)):
            val = np.random.uniform(0, 1)
            if val <= prob:
                self.cells[i][ab] = 1
            # else:
            #     self.cells[i][ab] = -1


        return

    def _assign_all_antibodies(self):
        """
        Assigns all the antibodies according to the probabilities in prob_arr
        """

        for i in range(len(self.ab_percentage_arr)):
            self._assign_ab(i, self.ab_percentage_arr[i])


    def get_marker_ratio(self, ab_list):
        """
        For all the cells, count the ones with all the markers in ab_seq and return their ratios
        :param ab_list: Is in format [A,B,C] : means antibodies A, B, and C are present
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

        self.marker_ratios[key] = ratio
        return ratio

    # def plot_patient(self):
    #
    #     x_data = []
    #     y_data = []
    #
    #     for i in range(self.ab_cnt):
    #
    #         x_data.append(i)
    #         y_data.append(self.get_marker_ratio([i]))
    #
    #     plt.plot(x_data, y_data)
    #
    #     plt.xlabel('Antibody index')
    #     plt.ylabel('Marker ratio')
    #     plt.grid(True)
    #
    #     plt.show()


# pt = Patient(100, 242, [1,2,3,4], 0.8, 0.1)
# pt.plot_patient()
# print pt.cells
# print pt.get_marker_ratio([1,2], [])

# for c in pt.cells:
#     print c
