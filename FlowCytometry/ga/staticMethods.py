
from itertools import combinations
import numpy as np


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
        key = ''
        if len(child) > 0:
            for i in range(len(child) - 1):
                key += str(int(child[i])) + '-'
            key += str(int(child[len(child) - 1]))

        return key

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
    def generate_cross_over_indices(ind_len):
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
    def get_unique_combinations(ab_arr):
        """
        Returns all the unique combinations of elements given in the antibody array, e.g. [a,b,c,d]
        :param ab_arr: Unsorted antibody array of 4
        :return: a list of unique 2^4 combinations of  antibodies
        """

        ab_arr = np.sort(ab_arr)

        zero_arr = np.full(len(ab_arr), 0, dtype=np.int)

        indices = np.concatenate((zero_arr, ab_arr))
        combs = combinations(indices, len(ab_arr))  # returns them sorted

        # remove duplicates from combinations
        comb_arr = [c for c in combs]
        comb_arr_unique = list(set(comb_arr))

        ab_list = []
        for comb in comb_arr_unique:
            present_ab = [ab for ab in comb if ab != 0]
            ab_list.append(present_ab)

        return ab_list

    @staticmethod
    def find_gaussian_intersection(mu_std_list):

        mu_int = mu_std_list[0][0]
        std_int = mu_std_list[0][1]

        i = 1
        while i < len(mu_std_list):
            mu2 = mu_std_list[i][0]
            std2 = mu_std_list[i][1]
            mu_int = (std2 * mu_int + std_int * mu2) / (std_int + std2)
            std_int = 1 / std_int + 1 / std2
            i += 2

        return [mu_int, std_int]


