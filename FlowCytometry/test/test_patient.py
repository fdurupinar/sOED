import unittest
from ga.patient import Patient

ab_cnt = 5
cell_cnt = 10
markers_list = [[1, 2, 3], [4]]
mu_list = [0.5, 0.1]
sigma_list = [0.1, 0.1]


class TestPatient(unittest.TestCase):

    def test_init(self):
        p = Patient("nc", cell_cnt, ab_cnt, markers_list, mu_list, sigma_list)
        self.assertEqual(len(p.cells), cell_cnt)
        self.assertEqual(p.ab_cnt, ab_cnt)
        self.assertEqual(p.markers_list, markers_list)
        self.assertEqual(p.mu_list, mu_list)
        self.assertEqual(p.sigma_list, sigma_list)


    def test_fill_in_cells(self):
        p = Patient("nc",cell_cnt, ab_cnt, markers_list, mu_list,  sigma_list)

        self.assertFalse(p.is_marker_arr[0])
        self.assertTrue(p.is_marker_arr[1])
        self.assertTrue(p.is_marker_arr[2])
        self.assertTrue(p.is_marker_arr[3])
        self.assertTrue(p.is_marker_arr[4])

    def test_get_marker_ratio(self):
        mu_list2 = [1, 0.01]
        sigma_list2 = [0.0001, 0.0001]

        p = Patient("nc",cell_cnt, ab_cnt, markers_list, mu_list2, sigma_list2)

        # print p.is_marker_arr
        # print p.cells

        ratio1 = p.get_marker_count([1, 2, 3], False)
        ratio2 = p.get_marker_count([4], False)

        self.assertEqual(ratio1, 1)
        self.assertEqual(ratio2, 0)


if __name__ == '__main__':
    unittest.main()