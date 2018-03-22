from patient import Patient
import numpy as np


class PatientFactory:

    def __init__(self, cnt, ab_cnt, markers_list,  cell_cnt, mu_list, sigma_list):
        """
        Create cnt patients of type <type>
        :param type: "c" or "nc"
        :param cnt: patient count to create
        """

        self.patients = []

        for i in range(cnt):
            p = Patient(cell_cnt, ab_cnt, markers_list, mu_list, sigma_list)
            self.patients.append(p)

    def get_marker_mean_and_variance(self, ab_list, to_be_recorded):
        """
        Among all the patients in this group, find the mean and variance values
        :param ab_list:
        :param to_be_recorded:
        :return:
        """

        perc_arr = []
        for p in self.patients:
            perc_arr.append(p.get_marker_ratio(ab_list, to_be_recorded))

        mu = np.mean(perc_arr)
        sigma = np.std(perc_arr)

        return [mu, sigma]


