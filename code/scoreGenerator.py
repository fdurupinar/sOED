import numpy as np
from itertools import combinations

#A,B,C,D combinations [-,-,-,-],[-,-,-,+], [-,-,+,-]. [-,-,+,+], [-,+,-,-]. [-,+,-,+]. [-,+,+,-], [-,+,+,+].
#[+,-,-,-],[+,-,-,+], [+,-,+,-]. [+,-,+,+], [+,+,-,-]. [+,+,-,+]. [+,+,+,-], [+,+,+,+]





class ScoreGenerator:

    def __init__(self, antibody_cnt, nc_cnt, c_cnt):

        self.patient_c = np.empty([nc_cnt, antibody_cnt, antibody_cnt, antibody_cnt, antibody_cnt])
        self.patient_c.fill(-1)
        self.patient_nc = np.empty([c_cnt, antibody_cnt, antibody_cnt, antibody_cnt, antibody_cnt])
        self.patient_nc.fill(-1)

        self.c_cnt = c_cnt
        self.nc_cnt = nc_cnt


    def draw_patient_percentage(self):
        """
        This will be  acquired by experiments
        :return: percentage of markers for 16 combinations
        """
        return np.random.uniform(0, 1)

    def assign_percentages_for_single_patient(self,  patient_arr, patient_ind, ab_combinations, is_c):
        """
        Computes the marker percentages for antibody combinations for a single patient
        :param patient_arr:
        :param patient_ind:
        :param ab_combinations: All combinations of the 4 antibodies
        :param is_c: boolean to return if patient is cancer
        """
        # TODO: percentages may be different whether patient is c or nc
        for c in ab_combinations:
            patient_arr[patient_ind][c[0]][c[1]][c[2]][c[3]] = self.draw_patient_percentage()

    def assign_percentages_for_all_patients(self, ab_arr):
        """
        Computes the marker percentages for antibody combinations for all cancer and non-cancer patients
        :param ab_arr: unsorted antibody array of 4 antibodies
        """

        ab_combinations = self.get_unique_combinations(ab_arr)

        for i in range(len(self.patient_nc)):
            self.assign_percentages_for_single_patient(self.patient_nc, i, ab_combinations, False)

        for i in range(len(self.patient_c)):
            self.assign_percentages_for_single_patient(self.patient_nc, i, ab_combinations, True)

    def get_unique_combinations(self,ab_arr):
        """Returns all the unique combinations of elements given in the antibody sequence, including their negatives"""

        ab_arr = np.sort(ab_arr)

        indices = np.array([0, 0, 0, 0, ab_arr[0], ab_arr[1], ab_arr[2], ab_arr[3]])
        combs = combinations(indices, 4)  # returns them sorted

        # remove duplicates from combinations
        comb_arr = [c for c in combs]
        comb_arr_unique = list(set(comb_arr))

        return comb_arr_unique

    def compute_precision_for_ab_sequence(self, ab_seq, threshold):
        """
        Compute the score for a single ab sequence
        :param ab_seq:
        :param threshold:
        :return:
        """

        nc_prec_cnt = 0
        for nc in self.patient_nc:
            # make sure that the array has a meaningful value
            if nc[ab_seq[0]][ab_seq[1]][ab_seq[2]][ab_seq[3]] < 0: # if no value is assigned yet
                prec_nc = self.predict_percentage_for_ab_sequence(nc, ab_seq)
            else:
                prec_nc = nc[ab_seq[0]][ab_seq[1]][ab_seq[2]][ab_seq[3]]

            if prec_nc < threshold:
                nc_prec_cnt += 1

        c_prec_cnt = 0
        for c in self.patient_c:
            if c[ab_seq[0]][ab_seq[1]][ab_seq[2]][ab_seq[3]] < 0: # if no value is assigned yet
                prec_c = self.predict_percentage_for_ab_sequence(c, ab_seq)
            else:
                prec_c = c[ab_seq[0]][ab_seq[1]][ab_seq[2]][ab_seq[3]]
            if  prec_c >= threshold:
                c_prec_cnt += 1

        return (float(c_prec_cnt)/ self.c_cnt) + (float(nc_prec_cnt)/self.nc_cnt)

    def compute_max_precision_for_ab_combination(self, ab_arr, threshold):
        """
        Compute all the scores for 16 combinations and find the maximum
        :param ab_arr: 4 antibodies listed in an unsorted way
        :param threshold:
        :return: precision of ab combination
        """

        ab_combinations = self.get_unique_combinations(ab_arr)

        max_prec = -1000
        for comb in ab_combinations:
            prec = self.compute_precision_for_ab_sequence(comb, threshold)
            if prec > max_prec:
                max_prec = prec


        return max_prec

    def get_independent_groups(self, ab_seq):
        """
        e.g. with input (a,b,c,d), will return [[a, [b,c]], [a,b,c,d], [[ab], [cd]]]
        :param ab_seq: sequence of antibodies
        :return:
        """

        if len(ab_seq) <= 1:
            return [ab_seq]
        if len(ab_seq) == 2:
            return [[[ab_seq[0]], [ab_seq[1]]], [[ab_seq[0],ab_seq[1]]]]

        ab_seq_smaller = ab_seq[1:]

        groups_smaller = self.get_independent_groups(ab_seq_smaller)

        groups = []
        for i in range(len(groups_smaller)):
            group = groups_smaller[i]
            # append the first element to the group on its own
            group_new = list(group)  # copy the original to keep it unchanged
            group_new = [[ab_seq[0]]] + group_new
            groups.append(group_new)
            for j in range(len(group)):
                group_new = list(group)  # copy the original to keep it unchanged
                el = group[j]
                el_new = [ab_seq[0]] + el # appended to the element array
                group_new[j] = el_new # replace the group here
                groups.append(group_new)

        return groups

    def predict_percentage_for_group_intersection(self, patient, groups):
        """
        :param self:
        :param patient: patient_arr[patient_ind]
        :param groups: a double array to show antibody sequences, e.g. [[a,b], [c]]
        :return: Independence assumption: precision([a,b]) * precision([c])
        """

        prec = 1
        for group in groups:
            # format group into a sequence
            ab_seq = np.zeros(4)
            for i in range(len(group)):
                ab_seq[i] = group[i]
            ab_seq = np.sort(ab_seq)

            prec *= patient[ab_seq[0]][ab_seq[1]][ab_seq[2]][ab_seq[3]]

        return prec

    def predict_percentage_for_ab_sequence(self, patient, ab_seq):
        """
        Computes all possible independent groups and finds their average precision
        E.g. for sequence (a,b,c,d) it looks at 14 groups of a*b*c*d, a*bcd, ab*cd, ac*bd, ..., etc.
        :param patient: patient_arr[patient_ind]
        :param ab_seq:
        :return: mean value of all the possible groups with already measured values
        """


        groups = self.get_independent_groups(ab_seq)

        prec = 0
        for g in groups:
            prec += self.predict_percentage_for_group_intersection(patient, g)

        prec /= len(groups)

        return prec


# sg = ScoreGenerator(20,20,20)
# sg.assign_percentages_for_all_patients([1,2,3,4])
# prec1 = sg.compute_precision_for_ab_sequence([0,0,1,2], 0.2)
# prec2 = sg.compute_precision_for_ab_sequence([0,0,3,4], 0.2)
# inters = sg.predict_precision_for_group_intersection([[1,2], [3,4]], 0.2)
# print prec1*prec2 == inters
# print sg.predict_precision_for_ab_sequence([1,2,3,4],0.2)
# print sg.compute_max_precision_for_ab_combination([1,2,3,4], 0.2)
# print sg.compute_max_precision_for_ab_combination([1,2,3,4], 0.2)

# print sg.get_independent_groups([1,2,3,4])


