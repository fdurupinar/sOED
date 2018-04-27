from unittest import TestCase
from util.staticMethods import StaticMethods
import numpy as np


class TestStaticMethods(TestCase):

    def test_get_ab_key(self):

        key = StaticMethods.get_ab_key(np.array([1, 2]))
        self.assertEqual(key, "1-2")

        key = StaticMethods.get_ab_key(np.array([]))
        self.assertEqual(key, "")

        key = StaticMethods.get_ab_key(np.array([1, 3, 4]))
        self.assertEqual(key, "1-3-4")

    def test_generate_cross_over_indices_2_point(self):

        inds = StaticMethods.generate_cross_over_indices_2_point(4)

        self.assertEqual(len(inds), 36)

    def test_generate_cross_over_indices_1_point(self):

        inds = StaticMethods.generate_cross_over_indices_1_point(4)

        self.assertEqual(len(inds), 16)

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

    def test_find_intersection_probability(self):

        val0 = StaticMethods.find_intersection_probability(70, 80, 40, 100)
        val1 = StaticMethods.find_intersection_probability(70, 80, 50, 100)
        val2 = StaticMethods.find_intersection_probability(70, 80, 70, 100)
        val3 = StaticMethods.find_intersection_probability(70, 80, 60, 100)

        self.assertEqual(val0, 0)
        self.assertGreater(val1, 0)
        self.assertEqual(val2, 1.0)
        self.assertGreater(val3, 0)
        self.assertLess(val3, 1.0)

