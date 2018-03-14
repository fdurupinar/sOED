from unittest import TestCase
from ga import GASolver

max_generations = 1
population_size = 10
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


            print ga.population[0][i]






    # def test_cross_over_generation(self):
    #     self.fail()
    #
    # def test_is_child_unique(self):
    #     self.fail()
    #
    # def test_cross_over(self):
    #     self.fail()
    #
    # def test_mutate(self):
    #     self.fail()
    #
    # def test_increment_current_generation(self):
    #     self.fail()
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
