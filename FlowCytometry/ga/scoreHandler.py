import numpy as np
from scipy import stats
from itertools import combinations

from patientFactory import PatientFactory
from staticMethods import StaticMethods


class ScoreHandler:

    def __init__(self, ab_cnt, nc_cnt, c_cnt, markers_list,  cell_cnt, mu_list, sigma_list):

        self.measured_list = []  # already measured, sorted ab sequences

        self.measured_dict = {}

        self.c_cnt = c_cnt
        self.nc_cnt = nc_cnt

        pf_c = PatientFactory("c", self.c_cnt, ab_cnt, markers_list,  cell_cnt, mu_list, sigma_list)

        self.patients_c = pf_c.patients
        pf_nc = PatientFactory("nc", self.nc_cnt, ab_cnt, markers_list, cell_cnt,  mu_list, sigma_list)
        self.patients_nc = pf_nc.patients

    def is_measured(self, ab_arr):
        """
        Check of ab_arr is measured already
        :param ab_arr: sorted antibodies np array. Can be any length.
        :return:
        """

        key = StaticMethods.get_ab_key(ab_arr)
        return key in self.measured_dict

    def _add_measured(self, ab_arr):
        """
        After sorting it, add ab_arr to the measured list if it doesn't exist already
        :param ab_arr: Can be of any length and unsorted
        :return:
        """

        ab_arr = np.sort(ab_arr)

        key = StaticMethods.get_ab_key(ab_arr)
        self.measured_dict[key] = 0

    def update_measured(self, ab_arr):
        """
        Add all the combinations of antibodies in ab_arr to the measured list in order to keep track
        :param ab_arr: np array of all 4 antibodies
        :return:
        """

        for i in range(1, len(ab_arr) + 1):
            combs = combinations(ab_arr, i)
            comb_arr = [c for c in combs]
            comb_arr_unique = list(set(comb_arr))

            for c in comb_arr_unique:
                self._add_measured(c)

    def _get_independent_groups(self, ab_list):
        """
        e.g. with input (a,b,c,d), will return [[a, [b,c]], [a,b,c,d], [[ab], [cd]]]
        :param ab_list: a list of antibodies
        :return:
        """

        if len(ab_list) <= 1:
            return [[ab_list]]
        if len(ab_list) == 2:
            return [[[ab_list[0]], [ab_list[1]]], [[ab_list[0], ab_list[1]]]]

        ab_seq_smaller = ab_list[1:]

        groups_smaller = self._get_independent_groups(ab_seq_smaller)

        groups = []
        for i in range(len(groups_smaller)):
            group = groups_smaller[i]
            # append the first element to the group on its own
            group_new = list(group)  # copy the original to keep it unchanged
            group_new = [[ab_list[0]]] + group_new
            groups.append(group_new)
            for j in range(len(group)):
                group_new = list(group)  # copy the original to keep it unchanged
                el = group[j]
                el_new = [ab_list[0]] + el # appended to the element array
                group_new[j] = el_new # replace the group here
                groups.append(group_new)


        return groups

    def _draw_intersection_ratio(self, patient, marker_cnt_arr):
        """
        Draw an intersection value from a gaussian distribution
        :param marker_ratios:
        :return:
        """

        pdf = StaticMethods.create_intersection_distribution(marker_cnt_arr, patient.cell_cnt)
        val = float(StaticMethods.draw_intersection_value(pdf)) / patient.cell_cnt
        return val


    def _predict_percentage_for_group_intersection(self, patient, group):
        """
        :param self:
        :param patient: patient_arr[patient_ind]
        :param group: a double array to show antibody sequences, e.g. [[a,b], [c]]
        :return: Independence assumption: precision([a,b]) * precision([c])
        """

        perc = 1

        combined_el = np.array([])
        for el in group:
            combined_el = np.concatenate((combined_el, el))

        combined_el = np.sort(combined_el)
        if self.is_measured(combined_el):
            return 0

        # if a combined value is known, skip this group
        # e.g. if ab is known, no need to predict a * b
        marker_cnt_arr = []
        for el in group:

            # even if one value in the group is unknown, return 0
            el = np.sort(el)
            if not self.is_measured(el):
                return 0

            # must be already measured or predicted
            marker_cnt_arr.append(patient.get_marker_count(el, True))  # no need to consider the absent ones

        if len(marker_cnt_arr) != 2:
            return 0
        #
        # print group
        # print marker_cnt_arr
        perc = self._draw_intersection_ratio(patient, marker_cnt_arr)

        return perc

    def _predict_percentage_for_ab_list(self, patient, groups):
        """
        Computes all possible independent groups and finds their average precision
        E.g. for sequence (a,b,c) it looks at 5 groups of a*b*c, a*cd, ab*c, bc*a, abc etc.
        Doesn't try to make predictions for unknown values. It uses already measured values
        :param patient:
        :param groups: a double array to show antibody sequences, e.g. [[a,b], [c]]
        :return: mean value of all the possible groups with already measured values
        """

        perc = 0
        cnt = 0
        for g in groups:
            perc += self._predict_percentage_for_group_intersection(patient, g)
            if perc > 0:
                cnt += 1  # ignore the 0 values -- unmeasured ones

        if cnt > 0:
            perc /= cnt
        else:
            perc = 0

        return perc

    def _compute_precision_for_ab_list(self, ab_list):
        """
        Compute the score for a single ab sequence
        :param ab_list: should be sorted
        :return:
        """
        if self.is_measured(ab_list):
            group1 = [nc.get_marker_count(ab_list, True) for nc in self.patients_nc]
            group2 = [c.get_marker_count(ab_list, True) for c in self.patients_c]

        else:
            # Create groups of combinations
            if len(ab_list) == 0:
                return 0

            groups = self._get_independent_groups(ab_list)

            group1 = [self._predict_percentage_for_ab_list(nc, groups) for nc in self.patients_nc]
            group2 = [self._predict_percentage_for_ab_list(c, groups) for c in self.patients_c]

        prec = abs(stats.ttest_ind(group1, group2)[0])
        # print group1
        # print group2
        # print prec

        if prec != prec:  # tests for NaN
            prec = 0
        return prec


    def compute_max_precision_for_ab_combination(self, ab_arr):
        """
        Compute all the scores for 16 combinations and find the maximum
        :param ab_arr: 4 antibodies listed in an unsorted way
        :return: precision of ab combination
        """

        ab_combinations = StaticMethods.get_unique_combinations(ab_arr)

        max_prec = -1000
        for comb in ab_combinations:
            comb = np.sort(comb)
            prec = self._compute_precision_for_ab_list(comb)
            if prec > max_prec:
                max_prec = prec
        return max_prec

    #########################################################################
    #  DEBUGGING METHODS
    ########################################################################

    def _compute_unmeasured_precision_for_ab_list(self, ab_list):
        """
        Compute the score for a single ab sequence
        :param ab_list: should be sorted
        :return:
        """

        if len(ab_list) == 0:
            return 0
        group1 = [nc.get_marker_count(ab_list, False) for nc in self.patients_nc]
        group2 = [c.get_marker_count(ab_list, False) for c in self.patients_c]
        prec = abs(stats.ttest_ind(group1, group2)[0])

        # print group1
        # print group2
        # print prec

        if prec != prec:  # NaN
            prec = 0

        return prec

    def compute_max_possible_precision(self, ab_arr):
        ab_combinations = StaticMethods.get_unique_combinations(ab_arr)

        print "Marker precisions"

        max_prec = -1000
        for comb in ab_combinations:
            comb = np.sort(comb)
            prec = self._compute_unmeasured_precision_for_ab_list(comb)

            str_comb = ', '.join(str(i) for i in comb)
            print '[' + str_comb + "] " + str(prec)
            if prec > max_prec:
                max_prec = prec

        return max_prec


# sg = ScoreHandler(20,20,20)
# print sg.get_independent_groups([1,2,3,4])
# print sg.get_unique_combinations([1,2,3,4])
