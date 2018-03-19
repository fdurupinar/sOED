import unittest
from patient import Patient
import numpy as np

ab_cnt = 5
cell_cnt = 10
markers = [1, 2]
mu = 0.5
sigma = 0.1


class TestPatient(unittest.TestCase):

    def test_init(self):
        p = Patient(cell_cnt, ab_cnt, markers, mu, sigma)
        self.assertEqual(len(p.cells), cell_cnt)
        self.assertEqual(p.ab_cnt, ab_cnt)
        self.assertEqual(p.markers, markers)
        self.assertEqual(p.mu, mu)
        self.assertEqual(p.sigma, sigma)

    def test_compute_ab_percentages(self):
        p = Patient(cell_cnt, ab_cnt, markers, mu, sigma)

        p._compute_ab_percentages()

        for i in range(ab_cnt):
            self.assertGreaterEqual(p.ab_percentage_arr[i], 0)
            self.assertLessEqual(p.ab_percentage_arr[i], 1)

    def test_assign_ab(self):
        p = Patient(cell_cnt, ab_cnt, markers, mu, sigma)

        #  make sure all the cells for antibody "1" are assigned to 1
        p._assign_ab(1, 1)

        for c in p.cells:
            self.assertEqual(c[1], 1)

    def test_assign_all_antibodies(self):
        p = Patient(cell_cnt, ab_cnt, markers, mu, sigma)
        p._assign_all_antibodies()

        for c in p.cells:
            for i in range(ab_cnt):
                self.assertGreaterEqual(c[i], 0)
                self.assertLessEqual(c[i], 1)

    def test_get_marker_ratio(self):
        mu = 1
        sigma = 0.0001
        p = Patient(cell_cnt, ab_cnt, markers, mu, sigma)

        ratio1 = p.get_marker_ratio([1], [])
        ratio2 = p.get_marker_ratio([2], [])
        ratio3 = p.get_marker_ratio([1, 2], [])
        ratio4 = p.get_marker_ratio([1], [2])

        self.assertEqual(ratio1, 1)
        self.assertEqual(ratio2, 1)
        self.assertEqual(ratio3, 1)
        self.assertEqual(ratio4, 0)


if __name__ == '__main__':
    unittest.main()