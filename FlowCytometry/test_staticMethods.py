from unittest import TestCase
from staticMethods import StaticMethods
import numpy as np

class TestStaticMethods(TestCase):

    def test_get_ab_key(self):

        key = StaticMethods.get_ab_key(np.array([1, 2]))
        self.assertEqual(key, "1-2")

        key = StaticMethods.get_ab_key(np.array([]))
        self.assertEqual(key, "")

        key = StaticMethods.get_ab_key(np.array([1, 3, 4]))
        self.assertEqual(key, "1-3-4")

    def test_generate_cross_over_indices(self):

        inds = StaticMethods.generate_cross_over_indices(4)

        self.assertEqual(len(inds), 36)

    def test_get_unique_combinations(self):

        ab_arr = np.array([4, 2, 1, 3])

        combs = StaticMethods.get_unique_combinations(ab_arr)

        # there are 16 combinations of an array of length 4
        self.assertEqual(len(combs), 16)

        self.assertTrue([1, 2] in combs)
        self.assertTrue([1, 2, 3, 4] in combs)
        self.assertTrue([] in combs)
        self.assertTrue([1, 3] in combs)
        self.assertTrue([1] in combs)
        self.assertTrue([2, 3, 4] in combs)

        self.assertFalse([2, 1] in combs, "should be sorted")
