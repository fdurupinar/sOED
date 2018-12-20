from unittest import TestCase
from ga.scoreHandlerHypothetical import ScoreHandlerHypothetical
import numpy as np

cell_cnt = 100

c_mu = 0.7
c_sigma = 0.1

nc_mu = 0.3
nc_sigma = 0.1

ab_cnt = 5
c_cnt = 10
nc_cnt = 10


markers_list = [[1, 2, 3, 4]]
mu_list = [0.99]
sigma_list = [0.1]


class TestScoreHandler(TestCase):

    def test_init(self):
        sh = ScoreHandlerHypothetical(ab_cnt, nc_cnt, c_cnt, markers_list, cell_cnt, mu_list,
                                      sigma_list)
        self.assertEqual(sh.c_cnt, c_cnt)
        self.assertEqual(sh.nc_cnt, nc_cnt)

        self.assertEqual(len(sh.patients_nc), nc_cnt)
        self.assertEqual(len(sh.patients_c), c_cnt)

    def test_is_measured(self):
        sh = ScoreHandlerHypothetical(ab_cnt, nc_cnt, c_cnt, markers_list, cell_cnt, mu_list,
                                      sigma_list)

        arr1 = np.array([1, 2])
        self.assertFalse(sh.is_measured(arr1))

        sh.measured_dict["1-2"] = 0
        self.assertTrue(sh.is_measured(arr1))

        arr2 = np.array([2, 1])
        self.assertFalse(sh.is_measured(arr2))

    def test__add_measured(self):
        sh = ScoreHandlerHypothetical(ab_cnt, nc_cnt, c_cnt, markers_list, cell_cnt, mu_list,
                                      sigma_list)

        arr1 = np.array([1, 2])
        sh._add_measured(arr1)

        arr2 = np.array([3, 4])
        sh._add_measured(arr2)

        self.assertTrue(sh.is_measured(arr1))
        self.assertTrue(sh.is_measured(arr2))

    def test_update_measured(self):
        sh = ScoreHandlerHypothetical(ab_cnt, nc_cnt, c_cnt, markers_list, cell_cnt, mu_list,
                                      sigma_list)

        ab_arr = np.array([4, 2, 1, 3])

        sh.update_measured(ab_arr)

        self.assertTrue(sh.is_measured(np.array([1, 2])))
        self.assertTrue(sh.is_measured(np.array([1, 2, 3, 4])))
        self.assertTrue(sh.is_measured(np.array([1, 3, 4])))
        self.assertTrue(sh.is_measured(np.array([2,  4])))
        self.assertTrue(sh.is_measured(np.array([3])))
        self.assertFalse(sh.is_measured(np.array([4, 3, 2, 1])), "should be sorted")

    def test_get_independent_groups(self):
        sh = ScoreHandlerHypothetical(ab_cnt, nc_cnt, c_cnt, markers_list, cell_cnt, mu_list,
                                      sigma_list)

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
        sh = ScoreHandlerHypothetical(10, nc_cnt, c_cnt, markers_list, cell_cnt, mu_list,
                                      sigma_list)

        sh._add_measured([5])
        sh._add_measured([6])
        sh._add_measured([7])
        sh._add_measured([8])

        patient = sh.patients_c[0]

        p1 = patient.get_marker_count([5], False)
        p2 = patient.get_marker_count([6], False)
        p3 = patient.get_marker_count([7], False)
        p4 = patient.get_marker_count([8], False)

        perc = sh._predict_percentage_for_group_intersection(patient, [[5], [6], [7], [8]])

        self.assertEqual(perc, p1 * p2 * p3 * p4)

        sh._add_measured([5, 6])
        sh._add_measured([7, 8])

        p12 = patient.get_marker_count([5, 6], False)
        p34 = patient.get_marker_count([7, 8], False)

        perc = sh._predict_percentage_for_group_intersection(patient, [[5, 6], [7, 8]])
        self.assertAlmostEqual(perc, p12 * p34)

    def test_predict_percentage_for_ab_list(self):
        sh = ScoreHandlerHypothetical(ab_cnt, nc_cnt, c_cnt, markers_list, cell_cnt, mu_list,
                                      sigma_list)

        # before updating measured values
        patient = sh.patients_c[0]
        perc1 = sh._predict_percentage_for_group_intersection(patient, [[1], [2]])
        perc2 = sh._predict_percentage_for_group_intersection(patient, [[1, 2]])
        perc = sh._compute_precision_for_ab_list([1, 2])

        self.assertEqual(perc, (perc1 + perc2)/2)

        sh.update_measured([1, 2])

        # after updating measured values
        perc1 = sh._predict_percentage_for_group_intersection(patient, [[1], [2]])
        # Because it doesn't cover maximally
        self.assertEqual(perc1, 0)

        perc2 = sh._predict_percentage_for_group_intersection(patient, [[1, 2]])
        # Because it is already measured
        self.assertEqual(perc2, 0)

    def test_compute_precision_for_ab_list(self):
        sh = ScoreHandlerHypothetical(10, 1, 1, markers_list, cell_cnt, mu_list,
                                      sigma_list)

        # already measured
        sh.update_measured([5, 6, 7, 8])

        prec = sh._compute_precision_for_ab_list([5, 6, 7, 8])
        self.assertEqual(prec, 0)

        prec1 = sh._compute_precision_for_ab_list([1, 2, 7, 8])
        prec2 = sh._compute_precision_for_ab_list([1, 2, 5, 6])
        self.assertEqual(prec1, prec2)

    def test_compute_max_precision_for_ab_combination(self):
        sh = ScoreHandlerHypothetical(5, 2, 2, markers_list, 10, [1],
                                      [0])

        # already measured
        sh.update_measured([1, 2, 3, 4])  # markers list

        prec1 = sh.compute_max_precision_for_ab_combination(np.array([1, 2, 3, 4]))
        prec2 = sh.compute_max_precision_for_ab_combination(np.array([1, 2, 4, 3]))
        self.assertEqual(prec1, prec2)

        prec = sh.compute_max_precision_for_ab_combination(np.array([1, 2, 3, 4]))

        self.assertGreaterEqual(prec, 1000) # inf

    def test__draw_ratio_intersection(self):
        sh = ScoreHandlerHypothetical(ab_cnt, nc_cnt, c_cnt, markers_list, cell_cnt, mu_list,
                                      sigma_list)

        marker_cnt_arr = [70, 80]

        val = sh._draw_intersection_ratio(sh.patients_c[0], marker_cnt_arr)

        self.assertGreaterEqual(val, 0.5)
        self.assertLessEqual(val, 0.7)

        marker_cnt_arr = [30, 40]

        val = sh._draw_intersection_ratio(sh.patients_c[0],marker_cnt_arr)

        self.assertGreaterEqual(val, 0)
        self.assertLessEqual(val, 0.3)
