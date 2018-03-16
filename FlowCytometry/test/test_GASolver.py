from unittest import TestCase
from ga import GASolver
import numpy as np


max_generations = 10
population_size = 20
antibody_cnt = 20
nc_cnt = 5
c_cnt = 5
markers = [1, 2, 3, 4]
c_threshold = 0.4

class TestGASolver(TestCase):

    def test_init(self):
        ga = GASolver(max_generations, population_size, antibody_cnt, nc_cnt, c_cnt, markers, c_threshold)

        self.assertEqual(ga.max_generations, max_generations)
        self.assertEqual(ga.population_size, population_size)
        self.assertEqual(ga.antibody_cnt, antibody_cnt)
        self.assertEqual(ga.nc_cnt, nc_cnt)
        self.assertEqual(ga.c_cnt, c_cnt)
        self.assertEqual(ga.markers, markers)
        self.assertEqual(ga.c_threshold, c_threshold)
        self.assertEqual(ga.total_population, population_size)

    def test_create_random_population(self):
        ga = GASolver(max_generations, population_size, antibody_cnt, nc_cnt, c_cnt, markers, c_threshold)
        ga.create_random_population()

        for i in range(population_size - 1):
            self.assertNotEqual(ga.population[0][i].tolist(), ga.population[0][i + 1].tolist(), "ab groups are different")
            for j in range(len(ga.population[0][i])-1):  # which is 4 -- 4 antibodies
                self.assertNotEqual(ga.population[0][i][j], ga.population[0][i][j+1], "ab's are different in each group")
                self.assertLess(ga.population[0][i][j], antibody_cnt, "antibody ids are correct")
                self.assertGreaterEqual(ga.population[0][i][j], 0, "antibody ids are correct")

    def test_get_fitness_key(self):
        ga = GASolver(max_generations, population_size, antibody_cnt, nc_cnt, c_cnt, markers, c_threshold)

        key = ga.get_fitness_key(np.array([1, 2, 3, 4]))
        self.assertEqual(key, "1-2-3-4")

        key = ga.get_fitness_key(np.array([0, 20, 15, 14]))
        self.assertEqual(key, "0-20-15-14")

        key = ga.get_fitness_key(np.array([0, 20, 15, 14, 3, 5, 7, 7]))
        self.assertEqual(key, "0-20-15-14-3-5-7-7")

    def test_get_fitness_value(self):
        ga = GASolver(max_generations, population_size, antibody_cnt, nc_cnt, c_cnt, markers, c_threshold)

        val = ga.get_fitness_value([1, 2, 3, 4])
        self.assertEqual(val, 0, "fitness not updated yet")

        ga.fitness["1-2-3-4"] = 0.5

        val = ga.get_fitness_value([1, 2, 3, 4])
        self.assertAlmostEqual(val, 0.5, "fitness updated")


    def test_is_child_diverse(self):

        ga = GASolver(max_generations, population_size, antibody_cnt, nc_cnt, c_cnt, markers, c_threshold)

        arr = np.array([1, 2, 3, 4])
        res = ga._is_child_diverse(arr)
        self.assertTrue(res)

        arr = np.array([1, 2, 3, 4, 1])
        res = ga._is_child_diverse(arr)
        self.assertFalse(res)

        arr = np.array([1, 2, 3, 4, 3, 1])
        res = ga._is_child_diverse(arr)
        self.assertFalse(res)

    def test_is_child_unique(self):
        ga = GASolver(max_generations, population_size, antibody_cnt, nc_cnt, c_cnt, markers, c_threshold)

        arr = np.array([1, 2, 3, 4])
        res = ga.is_child_unique(0, arr)
        self.assertTrue(res)

        arr = np.array([1, 2, 3, 4])
        ga.population[0][8] = arr

        res = ga.is_child_unique(0, arr)
        self.assertFalse(res)

    def test_generage_cross_over_indices(self):
        ga = GASolver(max_generations, population_size, antibody_cnt, nc_cnt, c_cnt, markers, c_threshold)

        inds = ga._generate_cross_over_indices(4)

        self.assertEqual(len(inds), 36)

    def test_cross_over(self):

        ga = GASolver(max_generations, population_size, antibody_cnt, nc_cnt, c_cnt, markers, c_threshold)

        group1 = np.array([1, 2, 3, 4])

        group2 = np.array([1, 6, 4, 8])

        child = ga.cross_over(group1, group2)


        self.assertTrue(child.tolist() != group1.tolist() and child.tolist() != group2.tolist())

        self.assertTrue(child[0] != child[1] and child[1] != child[2] and child[2] != child[3])

        self.assertTrue(child[0] < child[1] and child[1] < child[2] and child[2] < child[3])

    def test_mutate(self):

        ga = GASolver(max_generations, population_size, antibody_cnt, nc_cnt, c_cnt, markers, c_threshold)

        child = np.array([1, 2, 3, 4])
        mutant = ga.mutate(child)

        self.assertTrue(mutant[0] != mutant[1] and mutant[1] != mutant[2] and mutant[2] != mutant[3])
        self.assertTrue(mutant[0] < mutant[1] and mutant[1] < mutant[2] and mutant[2] < mutant[3])
        self.assertTrue(mutant[0] != child[0] or mutant[1] != child[1] or mutant[2] != child[2] or mutant[3] != child[3])

    def test_survive_n_fittest(self):

        ga = GASolver(max_generations, population_size, antibody_cnt, nc_cnt, c_cnt, markers, c_threshold)

        ga.create_random_population()
        fit_cnt = ga.survive_n_fittest(0, 2)

        self.assertNotEqual(ga.population[1][0].tolist(), [0, 0, 0, 0])
        self.assertNotEqual(ga.population[1][1].tolist(), [0, 0, 0, 0])
        self.assertEqual(ga.population[1][fit_cnt+1].tolist(), [0, 0, 0, 0])
        self.assertEqual(ga.population[1][fit_cnt+2].tolist(), [0, 0, 0, 0])

    def test_find_max_fitness_and_child(self):
        ga = GASolver(max_generations, population_size, antibody_cnt, nc_cnt, c_cnt, markers, c_threshold)
        ga.create_random_population()

        # don't include measured. All the population is measured so return -1
        val = ga.find_max_fitness_and_child(0, False)

        self.assertEqual(val['fitness'], -1)
        self.assertEqual(val['child'].tolist(), ga.population[0][0].tolist())

        # include measured
        val = ga.find_max_fitness_and_child(0, True)
        self.assertGreaterEqual(val['fitness'], 0)


        ga.cross_over_generation(0,ga.population_size/2)

        # don't include measured. Some are unmeasured so must be > 0
        val = ga.find_max_fitness_and_child(0, False)
        self.assertGreaterEqual(val['fitness'], 0)

        val = ga.find_max_fitness_and_child(0, True)
        self.assertGreaterEqual(val['fitness'], 0)

    def test_cross_over_generation(self):

        ga = GASolver(max_generations, population_size, antibody_cnt, nc_cnt, c_cnt, markers, c_threshold)


        ga.create_random_population()
        ga.survive_n_fittest(0, population_size / 2)

        start_ind = population_size / 2
        new_start_ind = ga.cross_over_generation(1, start_ind)
        self.assertEqual(ga.total_population, population_size * 5/4)
        self.assertEqual(new_start_ind, start_ind * 3/2)

        # half is zero
        self.assertTrue([0, 0, 0, 0] not in ga.population[0].tolist())
        x = [el for el in ga.population[1] if el.tolist() == [0, 0, 0, 0]]
        self.assertEqual(len(x), population_size/4)

    def test_mutate_generation(self):
        ga = GASolver(max_generations, population_size, antibody_cnt, nc_cnt, c_cnt, markers, c_threshold)

        np.random.rand(1)[0]
        ga.create_random_population()
        pop1 = np.copy(ga.population[0])

        ga.mutate_generation(0)
        pop2 = ga.population[1]

        self.assertNotEqual(pop1.tolist(), pop2.tolist())














        # def test_cross_over_generation(self):
    #     self.fail()
    #

    #

    #
    # def test_mutate(self):
    #     self.fail()
    #

    #
    # def test_mutate_generation(self):
    #     self.fail()
    #
    # def test_update_fitness(self):
    #     self.fail()
    #
    # def test_survive_n_fittest(self):
    #     self.fail()
    #
    # def test_run_simulation(self):
    #     self.fail()
