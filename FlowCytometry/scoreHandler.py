import numpy as np
from itertools import combinations

from patientFactory import PatientFactory


#A,B,C,D combinations [-,-,-,-],[-,-,-,+], [-,-,+,-]. [-,-,+,+], [-,+,-,-]. [-,+,-,+]. [-,+,+,-], [-,+,+,+].
#[+,-,-,-],[+,-,-,+], [+,-,+,-]. [+,-,+,+], [+,+,-,-]. [+,+,-,+]. [+,+,+,-], [+,+,+,+]


class ScoreHandler:

    def __init__(self, ab_cnt, nc_cnt, c_cnt, markers, cell_cnt, c_mu, c_sigma, nc_mu, nc_sigma):

        self.measured_list = []  # already measured, sorted ab sequences

        self.measured_dict = {}

        self.c_cnt = c_cnt
        self.nc_cnt = nc_cnt

        pf_c = PatientFactory("c", ab_cnt, markers, self.c_cnt, cell_cnt, c_mu, c_sigma, nc_mu, nc_sigma)
        self.patients_c = pf_c.patients
        pf_nc = PatientFactory("nc", ab_cnt, markers, self.nc_cnt, cell_cnt, c_mu, c_sigma, nc_mu, nc_sigma)
        self.patients_nc = pf_nc.patients


    @staticmethod
    def get_ab_key(child):
        """
        Encode child list as a string
        :param child: sorted np array of 4
        :return:

        """
        key = ''
        if len(child) > 0:
            for i in range(len(child)-1):
                key += str(int(child[i])) + '-'
            key += str(int(child[len(child)-1]))

        return key

    def is_measured(self, ab_arr):
        """
        Check of ab_arr is measured already
        :param ab_arr: sorted antibodies np array. Can be any length.
        :return:
        """

        key = self.get_ab_key(ab_arr)
        return key in self.measured_dict
        # return  ab_arr.tolist() in self.measured_list

    def _add_measured(self, ab_arr):
        """
        After sorting it, add ab_arr to the measured list if it doesn't exist already
        :param ab_arr:
        :return:
        """

        ab_arr = np.sort(ab_arr)

        key = self.get_ab_key(ab_arr)
        self.measured_dict[key] = 0

        # if not self.is_measured(ab_arr):
        #     self.measured_list.append(ab_arr.tolist())

    def _add_measurement(self, ab_arr, val):
        ab_arr = np.sort(ab_arr)

        key = self.get_ab_key(ab_arr)
        self.measured_dict[key] = val

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

    @staticmethod
    def _get_unique_combinations(ab_arr):
        """
        Returns all the unique combinations of elements given in the antibody array, e.g. [a,b,c,d]
        :param ab_arr: Unsorted antibody array of 4
        :return: a dict list of unique 2^4 combinations of present and absent antibodies [{'p':[a,b], 'a':[c,d]}, ...]
        """

        ab_arr = np.sort(ab_arr)

        zero_arr = np.zeros(len(ab_arr))

        indices = np.concatenate((zero_arr, ab_arr))
        combs = combinations(indices, len(ab_arr))  # returns them sorted

        # remove duplicates from combinations
        comb_arr = [c for c in combs]
        comb_arr_unique = list(set(comb_arr))

        ab_list = []
        for comb in comb_arr_unique:
            present_ab = [ab for ab in comb if ab != 0]
            absent_ab = [ab for ab in ab_arr if ab not in present_ab]
            ab_list.append({'p': present_ab, 'a': absent_ab})

        return ab_list

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
        # if combined_el.tolist() in self.measured_list: # this must have been observed before
        if self.is_measured(combined_el):
            return 0

        # if a combined value is known, skip this group
        # e.g. if ab is known, no need to predict a * b
        for el in group:

            # even if one value in the group is unknown, return 0
            el = np.sort(el)
            # if el.tolist() not in self.measured_list:
            if not self.is_measured(el):
                return 0

            #must be already measured or predicted
            perc *= patient.get_marker_ratio(el, [])  # no need to consider the absent ones


        return perc

    def _predict_percentage_for_ab_list(self, patient, ab_list):
        """
        Computes all possible independent groups and finds their average precision
        E.g. for sequence (a,b,c) it looks at 5 groups of a*b*c, a*cd, ab*c, bc*a, abc etc.
        Doesn't try to make predictions for unknown values. It uses already measured values
        :param patient:
        :param ab_list: could be any length
        :return: mean value of all the possible groups with already measured values
        """

        if len(ab_list) == 0:
            return 0

        # Create groups of combinations
        groups = self._get_independent_groups(ab_list)

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

    def _compute_precision_for_ab_list(self, present_ab_list, absent_ab_list, threshold):
        """
        Compute the score for a single ab sequence
        :param present_ab_list:
        :param absent_ab_list:
        :param threshold:
        :return:
        """
        nc_marker_cnt = 0
        c_marker_cnt = 0

        full_ab_list = np.concatenate((present_ab_list, absent_ab_list))
        full_ab_list = np.sort(full_ab_list)

        # if full_ab_list.tolist() in self.measured_list:  # read from the experiments
        if self.is_measured(full_ab_list):
            for nc in self.patients_nc:
                if nc.get_marker_ratio(present_ab_list, absent_ab_list) < threshold:
                    nc_marker_cnt += 1

            for c in self.patients_c:
                if c.get_marker_ratio(present_ab_list, absent_ab_list) >= threshold:
                    c_marker_cnt += 1
        else:
            for nc in self.patients_nc:
                if self._predict_percentage_for_ab_list(nc, present_ab_list) < threshold:
                    nc_marker_cnt += 1

            for c in self.patients_c:
                if self._predict_percentage_for_ab_list(c, present_ab_list) >= threshold:
                    c_marker_cnt += 1

        prec = float(c_marker_cnt + nc_marker_cnt) / (self.c_cnt + self.nc_cnt)

        return prec


    def compute_max_precision_for_ab_combination(self, ab_arr, threshold):
        """
        Compute all the scores for 16 combinations and find the maximum
        :param ab_arr: 4 antibodies listed in an unsorted way
        :param threshold:
        :return: precision of ab combination
        """

        ab_combinations = self._get_unique_combinations(ab_arr)

        max_prec = -1000
        for comb in ab_combinations:
            prec = self._compute_precision_for_ab_list(comb['p'], comb['a'], threshold)
            if prec > max_prec:
                max_prec = prec

        return max_prec

# sg = ScoreHandler(20,20,20)
# print sg.get_independent_groups([1,2,3,4])
# print sg.get_unique_combinations([1,2,3,4])
