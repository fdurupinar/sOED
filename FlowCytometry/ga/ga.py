
import matplotlib.pyplot as plt
import numpy as np
from scoreHandler import ScoreHandler
from util.staticMethods import StaticMethods
import operator

ANTIBODY_CNT = 64  # 16 # 64  # 240

POPULATION_SIZE = ANTIBODY_CNT / 4  # 40  # make this divisible by 4
MAX_GENERATIONS = 250 - POPULATION_SIZE
NC_CNT = 5  # non-cancer patients
C_CNT = 5  # cancer patients


CELL_CNT = 100


MARKERS_LIST = [[0, 5, 10, 15]]
MU_LIST = [0.99]
STD_DEV_LIST = [0.01]

MUTATION_PROBABILITY = 0.1  # 0.03
CROSS_OVER_2_RATIO = 0.5
VISUALIZE_POPULATION = False


class GASolver:

    def __init__(self, max_generations, population_size, antibody_cnt, nc_cnt, c_cnt, markers_list,  cell_cnt,
                 mu_list, sigma_list):

        self.max_generations = max_generations
        self.population_size = population_size
        self.antibody_cnt = antibody_cnt
        self.nc_cnt = nc_cnt
        self.c_cnt = c_cnt

        # fill with integers
        self.population = np.full((max_generations, population_size, 4), 0,  dtype=np.int)

        self.fitness = {}

        self.total_population = population_size

        self.score_handler = ScoreHandler(antibody_cnt, nc_cnt, c_cnt, markers_list, cell_cnt,
                                          mu_list, sigma_list)

        self.cross_over_indices_2_point = StaticMethods.generate_cross_over_indices_2_point(4)
        self.cross_over_indices_1_point = StaticMethods.generate_cross_over_indices_1_point(4)

        self.max_possible_fitness = self.score_handler.compute_max_possible_precision(MARKERS_LIST[0])

        self.plot_methods = PlotMethods(self)
        print "Maximum possible fitness value is:"
        print self.max_possible_fitness

    def create_population(self):
        """
        Create an initial population with all the possible antibodies
        :return:
        """

        i = 0
        while i < ANTIBODY_CNT:
            ab_arr = np.array([i, i+1, i+2, i+3])
            self.population[0][i/4] = ab_arr
            self.score_handler.update_measured(ab_arr)
            i += 4

    def assign_fitness(self, generation):
        """
        once the percentages are assigned, compute the precision scores deriving unknown values as multiplied
        probabilities of known values
        :return:
        """
        for i in range(self.population_size):
            self.update_fitness(self.population[generation][i])

    def create_random_population(self):
        """
        Initial creation of a population of size self.population_size
        :return:
        """

        # draw 4 random values from 242 antibodies n times
        # self.population = np.random.randint(ANTIBODY_CNT-1, size=(n, 4))

        for i in range(self.population_size):
            ab_arr = np.random.permutation(np.arange(0, self.antibody_cnt))[0:4]

            # sort
            ab_arr = np.sort(ab_arr)
            while ab_arr.tolist() in self.population[0].tolist():
                ab_arr = np.random.permutation(np.arange(0, self.antibody_cnt))[:4]
                ab_arr = np.sort(ab_arr)

            # convert all elements to integers first
            self.population[0][i] = ab_arr

            # initial assignment of percentages for all the patients for the ab_arr and its combinations
            self.score_handler.update_measured(ab_arr)

        # once the percentages are assigned, compute the precision scores of unknown values drawing from a distribution

        for i in range(self.population_size):
            self.update_fitness(self.population[0][i])

    def _is_child_diverse(self, child):
        """
        Are all the elements of the child are different
        :param child:
        :return:
        """
        uchild = np.unique(child)
        return len(uchild) == len(child)

    def is_child_unique(self, generation, child):
        """
        Check if there is already a child with the given sequence in that generation
        :param generation:
        :param child: np array
        :return:
        """
        return child.tolist() not in self.population[generation].tolist()

    def cross_over(self, group1, group2):
        """
        Perform 1 or 2-point cross-over between two groups of antibodies of size 4
        :param group1:
        :param group2:
        val = np.random.rand()

        :return:
        """
        child = np.full(len(group1), 0,  dtype=np.int)

        val = np.random.uniform(0, 1)
        if val < CROSS_OVER_2_RATIO:  # 2-point cross-over chances are low
            inds_order = np.arange(len(self.cross_over_indices_2_point))
            np.random.shuffle(inds_order)

            for ind in inds_order:

                inds1 = self.cross_over_indices_2_point[ind][0]
                inds2 = self.cross_over_indices_2_point[ind][1]

                # first half
                child[0] = group1[inds1[0]]
                child[1] = group1[inds1[1]]

                # second half
                child[2] = group2[inds2[0]]
                child[3] = group2[inds2[1]]

                if self._is_child_diverse(child):
                    break
        else:
            inds_order = np.arange(len(self.cross_over_indices_1_point))
            np.random.shuffle(inds_order)

            child = np.array(group1, copy=True)
            for ind in inds_order:

                inds1 = self.cross_over_indices_2_point[ind][0]
                inds2 = self.cross_over_indices_2_point[ind][1]

                # first half
                child[inds1] = group2[inds2]

                if self._is_child_diverse(child):
                    break


        child = np.sort(child)

        return child

    def cross_over_generation(self, generation, end_ind, fill_ind):
        """
        Cross over gene groups with the given ratio
        :param generation: generation number to reproduce
        :param end_ind: last index to cross over
        :param fill_ind: starting index to fill in the new generation
        :return:
        """
        # randomly divide these into two
        co_inds1 = np.random.permutation(end_ind)[0:end_ind / 2]
        co_inds2 = np.setdiff1d(np.arange(end_ind), co_inds1)
        #
        # print "cross over indices"
        # print co_inds1
        # print co_inds2

        for i in range(len(co_inds1)):
            child = self.cross_over(self.population[generation][co_inds1[i]], self.population[generation][co_inds2[i]])

            # if cross-over fails to generate a proper child, generate a random child
            while not (self.is_child_unique(generation, child) and self._is_child_diverse(child)):
                child = np.random.permutation(self.antibody_cnt)[0:4]
                child = np.sort(child)

            self.population[generation][fill_ind] = child
            self.update_fitness(child)

            fill_ind += 1
            self.total_population += 1  # increment total population

        return fill_ind

    def mutate(self, child):
        """
        Randomly mutate one antibody in a group of 4
        :param child:
        :return:
        """

        mut_ind = np.random.randint(len(child))

        ab = np.random.randint(self.antibody_cnt)

        # keep on random selection if the current antibody already exists
        while ab in child:
            ab = np.random.randint(self.antibody_cnt)

        mutant = np.copy(child)
        mutant[mut_ind] = ab
        mutant = np.sort(mutant)

        return mutant

    def mutate_generation(self, generation, start_ind, end_ind):
        """Mutate some portion of the new generation"""

        for i in range(start_ind, end_ind):
            child = self.population[generation][i]

            mutation_chance = np.random.uniform(0, 1)
            if mutation_chance < MUTATION_PROBABILITY:
                mutated_child = self.mutate(child)

                # if child exists, generate another child
                while not self.is_child_unique(generation, mutated_child):
                    mutated_child = self.mutate(mutated_child)

                self.population[generation][i] = mutated_child


    def get_fitness_value(self, child):
        """
        Return fitness value for child
        :param child: sorted np array of 4
        :return:
        """
        val = 0

        if StaticMethods.get_ab_key(child) in self.fitness:
            val = self.fitness[StaticMethods.get_ab_key(child)]

        return val

    def update_fitness(self, child):
        """
        Compute the fitness of antibody quadruple based on existing values
        :param child:
        :return:
        """

        score = self.score_handler.compute_max_precision_for_ab_combination(child)

        # print str(child) + " " + str(score)
        key = StaticMethods.get_ab_key(child)

        self.fitness[key] = score

    def survive_n_fittest(self, generation, n):
        """Find and survive the n fittest combinations in the population-- will overwrite the rest"""

        # fitness = np.zeros(len(self.population[generation]))
        sorted_fitness = sorted(self.fitness.items(), key=operator.itemgetter(1))
        len_fitness = len(sorted_fitness) - 1

        for i in range(n):
            self.population[generation + 1][i] = StaticMethods.key_to_ab(sorted_fitness[len_fitness - i][0])

    def find_max_fitness_and_child(self, generation, include_measured):
        """
        Find unmeasured fittest child so that we can measure it
        :param include_measured:
        :param generation:
        :return:
        """

        max_fitness = -1
        fittest_child = self.population[generation][0]
        for child in self.population[generation]:
            fitness = self.get_fitness_value(child)
            if fitness > max_fitness and (include_measured or not self.score_handler.is_measured(child)):
                max_fitness = fitness
                fittest_child = child

        return {'child': fittest_child, 'fitness': max_fitness}

    def run_simulation(self, max_gen_cnt):
        """
        Run evolution for max_gen_cnt iterations
        :param max_gen_cnt:
        :return:
        """
        if VISUALIZE_POPULATION:
            plt.ion()

        # self.create_random_population()
        self.create_population()

        for gen in range(max_gen_cnt-1):
            if VISUALIZE_POPULATION:
                self.plot_methods.plot_generation(gen)
                # self.plot_methods.visualize_generation(gen)

            # updates the fitness values of the whole population
            self.assign_fitness(gen)

            # survive half of the generation and pass them to the next gen
            self.survive_n_fittest(gen, self.population_size / 2)

            # cross over the next generation half

            new_start_ind = self.cross_over_generation(gen + 1, self.population_size / 2, self.population_size / 2)
            # cross over twice to keep population constant
            self.cross_over_generation(gen + 1, self.population_size / 2,  new_start_ind)

            # mutate the new ones only, not the old ones
            self.mutate_generation(gen + 1, self.population_size / 2, self.population_size)

            # find fittest unmeasured child
            unmeasured_max_fitness = self.find_max_fitness_and_child(gen + 1, False)
            # find the fittest child
            total_max_fitness = self.find_max_fitness_and_child(gen + 1, True)

            if unmeasured_max_fitness['fitness'] < 0:
                print "no unmeasured fitness value found"
            else:
                print "unmeasured \t total"
                print str(gen) + "\t" + str(unmeasured_max_fitness) + "\t" + str(total_max_fitness)

            if total_max_fitness['fitness'] >= self.max_possible_fitness:
                print "Success: "
                print total_max_fitness['child']
                break
            else:  # update values for fittest unmeasured child through the experiment
                self.score_handler.update_measured(unmeasured_max_fitness['child'])


    def _print_generation(self, generation):
        for p in self.population[generation]:
            print {'fitness': self.get_fitness_value(p), 'child': p}


gs = GASolver(MAX_GENERATIONS, POPULATION_SIZE, ANTIBODY_CNT, NC_CNT, C_CNT, MARKERS_LIST, CELL_CNT, MU_LIST, STD_DEV_LIST)
gs.run_simulation(MAX_GENERATIONS)

