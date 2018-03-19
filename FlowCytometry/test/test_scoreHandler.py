from unittest import TestCase
from scoreHandler import ScoreHandler
from patient import Patient
import numpy as np

cell_cnt = 100

c_mu = 0.7
c_sigma = 0.1

nc_mu = 0.3
nc_sigma = 0.1

ab_cnt = 5
c_cnt = 10
nc_cnt = 10
markers = [1, 2, 3, 4]

class TestScoreHandler(TestCase):

    def test_init(self):
        sh = ScoreHandler(ab_cnt, c_cnt, nc_cnt, markers,  cell_cnt, c_mu, c_sigma, nc_mu, nc_sigma)
        self.assertEqual(sh.c_cnt, c_cnt)
        self.assertEqual(sh.nc_cnt, nc_cnt)

        self.assertEqual(len(sh.patients_nc), nc_cnt)
        self.assertEqual(len(sh.patients_c), c_cnt)

    def test_is_measured(self):
        sh = ScoreHandler(ab_cnt, c_cnt, nc_cnt, markers, cell_cnt, c_mu, c_sigma, nc_mu, nc_sigma)

        arr1 = np.array([1,2])
        self.assertFalse(sh.is_measured(arr1))

        sh.measured_list.append([1, 2])
        self.assertTrue(sh.is_measured(arr1))

        arr2 = np.array([2, 1])
        self.assertFalse(sh.is_measured(arr2))

    def test__add_measured(self):
        sh = ScoreHandler(ab_cnt, c_cnt, nc_cnt, markers, cell_cnt, c_mu, c_sigma, nc_mu, nc_sigma)

        arr1 = np.array([1, 2])
        sh._add_measured(arr1)

        arr2 = np.array([3, 4])
        sh._add_measured(arr2)

        self.assertTrue(sh.is_measured(arr1))
        self.assertTrue(sh.is_measured(arr2))

    def test_update_measured(self):
        sh = ScoreHandler(ab_cnt, c_cnt, nc_cnt, markers, cell_cnt, c_mu, c_sigma, nc_mu, nc_sigma)

        ab_arr = np.array([4,2, 1, 3])

        sh.update_measured(ab_arr)

        self.assertTrue(sh.is_measured(np.array([1, 2])))
        self.assertTrue(sh.is_measured(np.array([1, 2, 3, 4])))
        self.assertTrue(sh.is_measured(np.array([1, 3, 4])))
        self.assertTrue(sh.is_measured(np.array([2,  4])))
        self.assertTrue(sh.is_measured(np.array([3])))
        self.assertFalse(sh.is_measured(np.array([4, 3, 2, 1])), "should be sorted")

    def test_get_independent_groups(self):
        sh = ScoreHandler(ab_cnt, c_cnt, nc_cnt, markers, cell_cnt, c_mu, c_sigma, nc_mu, nc_sigma)

        groups = sh._get_independent_groups([1, 2])

        # there are 2 groups for an array of length 2
        self.assertEqual(len(groups), 2)
        self.assertIn([[1], [2]], groups)
        self.assertIn([[1, 2]], groups)
        self.assertNotIn([[1, 2, 3]], groups, "extra element")
        self.assertNotIn([[2, 1]], groups, "should be sorted")

        groups = sh._get_independent_groups([1, 2, 3, 4])

        # there are 15 groups for an array of length 4
        self.assertEqual(len(groups), 15)
        self.assertIn([[1], [2], [3, 4]], groups)
        self.assertIn([[1, 2], [3, 4]], groups)
        self.assertIn([[2, 3], [1, 4]], groups)
        self.assertIn([[1, 2, 3], [4]], groups)
        self.assertIn([[1, 2, 3, 4]], groups)
        self.assertIn([[1], [2], [3], [4]], groups)
        self.assertNotIn([[1, 2], [4]], groups, "missing element")
        self.assertNotIn([[1, 2], [1, 3, 4]], groups, "cannot include same element in different elements")
        self.assertNotIn([[4, 3, 2, 1]], groups, "should be sorted")
        self.assertNotIn([[4], [3], [2], [1]], groups, "should be sorted")

    def test_predict_percentage_for_group_intesection(self):
        sh = ScoreHandler(ab_cnt, c_cnt, nc_cnt, markers, cell_cnt, c_mu, c_sigma, nc_mu, nc_sigma)

        sh.update_measured([1, 2, 3, 4])
        patient = sh.patients_c[0]

        p1 = patient.get_marker_ratio([1], [])
        p2 = patient.get_marker_ratio([2], [])
        p3 = patient.get_marker_ratio([3], [])
        p4 = patient.get_marker_ratio([4], [])

        perc = sh._predict_percentage_for_group_intersection(patient, [[1], [2], [3], [4]])

        self.assertAlmostEqual(perc, p1 * p2 * p3 * p4)

        p12 = patient.get_marker_ratio([1, 2], [])
        p34 = patient.get_marker_ratio([3, 4], [])

        perc = sh._predict_percentage_for_group_intersection(patient, [[1, 2], [3, 4]])
        self.assertAlmostEqual(perc, p12 * p34)

        p1234 = patient.get_marker_ratio([1, 2, 3, 4], [])

        perc = sh._predict_percentage_for_group_intersection(patient, [[1, 2, 3, 4]])
        self.assertAlmostEqual(perc, p1234)

    def test_predict_percentage_for_ab_list(self):
        sh = ScoreHandler(ab_cnt, c_cnt, nc_cnt, [1, 2],  cell_cnt, c_mu, c_sigma, nc_mu, nc_sigma)

        sh.update_measured([1, 2])
        patient = sh.patients_c[0]

        perc1 = sh._predict_percentage_for_group_intersection(patient, [[1], [2]])
        perc2 = sh._predict_percentage_for_group_intersection(patient, [[1, 2]])
        perc = sh._predict_percentage_for_ab_list(patient, [1, 2])

        self.assertEqual(perc, (perc1 + perc2)/2)

    def test_compute_precision_for_ab_list(self):
        sh = ScoreHandler(ab_cnt, 1, 1, [1, 2],  cell_cnt, c_mu, c_sigma, nc_mu, nc_sigma)  # 1 patient each

        # already measured
        sh.update_measured([1, 2])

        prec = sh._compute_precision_for_ab_list([1, 2], [], 1)
        self.assertEqual(prec, 0.5)

        prec = sh._compute_precision_for_ab_list([1, 2], [], 0)
        self.assertEqual(prec, 0.5)

        prec = sh._compute_precision_for_ab_list([1, 2], [], 0.4)
        self.assertGreaterEqual(prec, 0.5)

    def test_compute_max_precision_for_ab_combination(self):
        sh = ScoreHandler(5, 1, 1, [1, 2, 3, 4],  cell_cnt, c_mu, c_sigma, nc_mu, nc_sigma)  # 1 patient each

        # already measured
        sh.update_measured([1, 2, 3, 4])

        prec = sh.compute_max_precision_for_ab_combination(np.array([1, 2, 4, 3]), 0)
        self.assertEqual(prec, 0.5)

        prec = sh.compute_max_precision_for_ab_combination(np.array([1, 2, 4, 3]), 1)
        self.assertEqual(prec, 0.5)

        prec = sh.compute_max_precision_for_ab_combination(np.array([1, 2, 4, 3]), 0.4)
        self.assertGreaterEqual(prec, 0.5)

        #  clear sh and customly assign patients
        sh = ScoreHandler(5, 2, 2, [1, 2],  cell_cnt, c_mu, c_sigma, nc_mu, nc_sigma)  # 2 patients each
        sh.update_measured([1, 2])
        #
        sh.patients_c = []
        sh.patients_nc = []
        sh.patients_nc.append(Patient(5, 10, [1, 2], 0, 0.001))  # non-cancer patients don't have these markers
        sh.patients_c.append(Patient(5, 10, [1, 2], 1, 0.001))  # cancer patients have these markers
        sh.patients_nc.append(Patient(5, 10, [1, 2], 0, 0.001))  # non-cancer patients don't have these markers
        sh.patients_c.append(Patient(5, 10, [1, 2], 1, 0.001))  # cancer patients have these markers

        prec = sh.compute_max_precision_for_ab_combination(np.array([1, 2]), 0.4)
        self.assertEqual(prec, 1)

        # test where 1 and 2 are all 1s for all patients, and others are all 0s for all patients

        sh.update_measured([1, 2, 3, 4])

        # make other cells 0
        for j in range(2):
            for i in range(5):
                sh.patients_nc[j].cells[i][0] = 0
                sh.patients_nc[j].cells[i][3] = 0
                sh.patients_nc[j].cells[i][4] = 0
                sh.patients_c[j].cells[i][0] = 0
                sh.patients_c[j].cells[i][3] = 0
                sh.patients_c[j].cells[i][4] = 0

                sh.patients_nc[j].cells[i][1] = 1
                sh.patients_nc[j].cells[i][2] = 1
                sh.patients_c[j].cells[i][1] = 1
                sh.patients_c[j].cells[i][2] = 1

        prec = sh.compute_max_precision_for_ab_combination(np.array([1, 2, 3, 4]), 0.4)

        self.assertEqual(prec, 0.5)

        # test 3/4 precision case
        for j in range(2):
            for i in range(5):
                sh.patients_nc[j].cells[i][0] = 0
                sh.patients_nc[j].cells[i][3] = 0
                sh.patients_nc[j].cells[i][4] = 0
                sh.patients_c[j].cells[i][0] = 0
                sh.patients_c[j].cells[i][3] = 0
                sh.patients_c[j].cells[i][4] = 0
                sh.patients_c[j].cells[i][1] = 1
                sh.patients_c[j].cells[i][2] = 1
                sh.patients_nc[j].cells[i][2] = 1

        for i in range(5):  # only 1 nc is correct
            sh.patients_nc[0].cells[i][1] = 0
            sh.patients_nc[1].cells[i][1] = 1

        prec = sh.compute_max_precision_for_ab_combination(np.array([1, 2, 3, 4]), 0.4)

        self.assertEqual(prec, 0.75)

