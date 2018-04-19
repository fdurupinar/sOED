
from itertools import combinations
import numpy as np
import scipy.stats as sp


class StaticMethods:

    def __init__(self):
        return

    @staticmethod
    def get_ab_key(child):
        """
        Encode child list as a string
        :param child: sorted np array of 4
        :return:

        """

        return '-'.join(str(c) for c in child)
        # key = ''
        # if len(child) > 0:
        #     for i in range(len(child) - 1):
        #         key += str(int(child[i])) + '-'
        #     key += str(int(child[len(child) - 1]))

        # return key

    @staticmethod
    def key_to_ab(key):
        """
        Decode key string into antibody indices
        :param key: str
        :return:

        """

        child = key.split("-")

        child = [int(c) for c in child]

        return child

    @staticmethod
    def generate_cross_over_indices_2_point(ind_len):
        """
        Generate an array of arrays for indices to cross-over
        :param ind_len:
        :return:
        """
        inds_group = []
        for i in range(0, ind_len):
            for j in range(i + 1, ind_len):
                for k in range(0, ind_len):
                    for l in range(k + 1, ind_len):
                        inds = [[i, j], [k, l]]
                        inds_group.append(inds)

        return inds_group

    @staticmethod
    def generate_cross_over_indices_1_point(ind_len):
        """
        Generate an array of arrays for indices to cross-over
        Generates all possible combinations to be randomly selected from
        :param ind_len:
        :return:
        """
        inds_group = []
        for i in range(0, ind_len):
            for k in range(0, ind_len):
                inds = [[i, k]]
                inds_group.append(inds)

        return inds_group

    @staticmethod
    def get_unique_combinations(ab_arr):
        """
        Returns all the unique combinations of elements given in the antibody array, e.g. [a,b,c,d]
        :param ab_arr: Unsorted antibody array of 4
        :return: a list of unique 2^4 combinations of  antibodies
        """

        ab_arr = np.sort(ab_arr)


        zero_arr = np.full(len(ab_arr), -1, dtype=np.int)

        indices = np.concatenate((zero_arr, ab_arr))
        combs = combinations(indices, len(ab_arr))  # returns them sorted

        # remove duplicates from combinations
        comb_arr = [c for c in combs]
        comb_arr_unique = list(set(comb_arr))

        ab_list = []
        for comb in comb_arr_unique:
            present_ab = [ab for ab in comb if ab != -1]
            ab_list.append(present_ab)

        return ab_list

    @staticmethod
    def create_intersection_distribution(marker_cnt_arr, cell_cnt):
        """
        Create a distribution with cnt1/cell_cnt, cnt2/cell_cnt ratio groups with intersection probabilities
        :param marker_cnt_arr:
        :param cell_cnt:
        :return: dist
        """
        cnt1 = marker_cnt_arr[0]
        cnt2 = marker_cnt_arr[1]

        start_val = np.max([0, cnt1 + cnt2 - cell_cnt])
        end_val = np.min([cnt1, cnt2])

        # cumulative distribution function
        cdf = []
        val = start_val
        x_range = []
        while val <= end_val:
            prob = StaticMethods.find_intersection_probability(cnt1, cnt2, val, cell_cnt)
            # print str(cnt1) + " " + str(cnt2) +  " " + str(val)  + " "  + str(cell_cnt) + " " + str(prob)
            cdf.append(prob)
            x_range.append(val)

            val += 1.0


        # probability distribution function
        # pdf = [0]
        # for i in range(len(cdf) - 1):
        #     pdf.append(cdf[i+1] - cdf[i])


        return {'x': x_range, 'y': cdf}

    @staticmethod
    def find_intersection_probability(cnt1, cnt2, cnt12, cell_cnt):
        """
        Find the probability that intersection of two groups is less than int12/cellCnt
        :param cnt1: Number of markers in group 1
        :param cnt2: Number of markers in group s
        :param cnt12: Number of markers in intersection
        :param cell_cnt: Total number of cells
        :return:
        """

        try:
            prob = sp.fisher_exact([[cnt12, cnt2 - cnt12], [cnt1 - cnt12, cell_cnt - (cnt1 + cnt2 - cnt12)]], 'less')
            return prob[1]
        except:
            return 0

    @staticmethod
    def draw_intersection_value(cdf):
        """
        Given a cumulative distribution function, draw a value randomly
        :param cdf: dictionary of {'x':[], 'y':[]} y's are the pdf values, x's are the overlap values to be returned
        :return:
        """

        y_vals = cdf['y']

        val = np.random.rand()

        for i in range(len(y_vals)-1):
            if val >= y_vals[i] and val <= y_vals[i+1]:
                return cdf['x'][i]

        return 0


# print StaticMethods.find_intersection_probability(70, 80, 50,  100)
#
#
# data = StaticMethods.create_intersection_distribution([70, 80], 100)
# print StaticMethods.draw_intersection_value(data)
# # #
# print data
# #
# import matplotlib.pyplot as plt
# plt.bar(data['x'], data['y'])
# plt.show()