
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
import matplotlib.animation as anim
from scoreHandler import ScoreHandler


ANTIBODY_CNT = 50
MUTATION_PROBABILITY = 0.03
POPULATION_SIZE = 40  # make this divisible by 4
MAX_GENERATIONS = 100
NC_CNT = 20  # non-cancer patients
C_CNT = 20  # cancer patients

C_THRESHOLD = 0.4


CELL_CNT = 100

CANCER_MU = 0.7
CANCER_STD_DEV = 0.1

NON_CANCER_MU = 0.3
NON_CANCER_STD_DEV = 0.1


MARKERS = [1, 5, 7, 8]

class GASolver:

    def __init__(self, max_generations, population_size, antibody_cnt, nc_cnt, c_cnt, markers, c_threshold, cell_cnt,
                 c_mu, c_sigma, nc_mu, nc_sigma):

        self.experiment_cnt = 0
        self.max_generations = max_generations
        self.population_size = population_size
        self.antibody_cnt = antibody_cnt
        self.nc_cnt = nc_cnt
        self.c_cnt = c_cnt
        self.markers = markers
        self.c_threshold = c_threshold

        # fill with integers
        self.population = np.full((max_generations, population_size, 4), 0,  dtype=np.int)

        self.fitness = {}

        self.total_population = population_size

        self.score_handler = ScoreHandler(antibody_cnt, nc_cnt, c_cnt, markers, cell_cnt, c_mu, c_sigma, nc_mu, nc_sigma)

        self.cross_over_indices = self._generate_cross_over_indices(4)

    @staticmethod
    def _generate_cross_over_indices(ind_len):
        """
        Generate an array of arrays for indices to cross-over
        :param ind_len:
        :return:
        """
        inds_group = []
        for i in range(0, ind_len):
            for j in range(i + 1, ind_len):
                for k in range(0, ind_len):
                    for l in range(k + 1, ind_len):
                        inds = [[i, j], [k, l]]
                        inds_group.append(inds)

        return inds_group


    def create_random_population(self):
        """
        Initial creation of a population of size self.population_size
        :return:
        """

        # draw 4 random values from 242 antibodies n times
        # self.population = np.random.randint(ANTIBODY_CNT-1, size=(n, 4))


        for i in range(self.population_size):
            ab_arr = np.random.permutation(np.arange(0, self.antibody_cnt))[0:4]

            # creating a random population without repeating elements is costly
            # sort
            ab_arr = np.sort(ab_arr)
            while ab_arr.tolist() in self.population[0].tolist():
                ab_arr = np.random.permutation(np.arange(0, self.antibody_cnt))[0:4]
                ab_arr = np.sort(ab_arr)

            # convert all elements to integer first
            self.population[0][i] = ab_arr

            # initial assignment of percentages for all the patients for the ab_arr and its combinations
            self.score_handler.update_measured(ab_arr)


        # once the percentages are assigned, compute the precision scores deriving unknown values as multiplied
        # probabilities of known values

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
        if child.tolist() in self.population[generation].tolist():
                return False
        else:
            return True



    def cross_over(self, group1, group2):
        """
        Perform 2-point cross-over between two groups of antibodies of size 4
        :param group1:
        :param group2:
        :return:
        """
        child = np.full(len(group1), 0,  dtype=np.int)

        inds_order = np.arange(len(self.cross_over_indices))
        np.random.shuffle(inds_order)

        for ind in inds_order:

            inds1 = self.cross_over_indices[ind][0]
            inds2 = self.cross_over_indices[ind][1]

            # first half
            child[0] = group1[inds1[0]]
            child[1] = group1[inds1[1]]

            # second half
            child[2] = group2[inds2[0]]
            child[3] = group2[inds2[1]]

            if self._is_child_diverse(child):
                break

        child = np.sort(child)

        return child

    def cross_over_generation(self, generation, start_ind):
        """
        Cross over gene groups with the given ratio
        :param generation: generation number to reproduce
        :param start_ind: starting index to fill in the new generation
        :return:
        """

        # randomly divide these into two
        co_inds1 = np.random.permutation(start_ind)[0:self.population_size / 4]
        co_inds2 = np.setdiff1d(np.arange(start_ind), co_inds1)

        for i in range(len(co_inds1)):

            child = self.cross_over(self.population[generation][co_inds1[i]], self.population[generation][co_inds2[i]])

            # if cross-over fails to generate a proper child, generate a random child
            if not self.is_child_unique(generation, child) or child.tolist() == [0, 0, 0, 0]:
                child = np.random.permutation(self.antibody_cnt)[0:4]
                while child.tolist() in self.population[generation].tolist():
                    child = np.random.permutation(self.antibody_cnt)[0:4]
                    child = np.sort(child)

            self.population[generation][start_ind] = child
            self.update_fitness(child)

            start_ind += 1
            self.total_population += 1  # increment total population

        return start_ind

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

        if self.score_handler.get_ab_key(child) in self.fitness:
            val = self.fitness[self.score_handler.get_ab_key(child)]

        return val

    def update_fitness(self, child):
        """
        Compute the fitness of antibody quadruple based on existing values
        :param child:
        :return:
        """

        score = self.score_handler.compute_max_precision_for_ab_combination(child, self.c_threshold)

        key = self.score_handler.get_ab_key(child)
        self.fitness[key] = score

    def survive_n_fittest(self, generation, n):
        """Find and survive the n fittest combinations in the population-- will overwrite the rest"""

        fitness = np.zeros(len(self.population[generation]))
        for i in range(len(self.population[generation])):
            child = self.population[generation][i]

            fitness[i] = self.get_fitness_value(child)  # self.compute_fitness(child)

        fitness = np.sort(fitness)
        threshold = fitness[len(fitness) - n]

        j = 0
        for i in range(len(self.population[generation])):
            child = self.population[generation][i]
            if self.get_fitness_value(child) >= threshold:
                self.population[generation + 1][j] = child
                j += 1

        return j

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

        self.create_random_population()

        for i in range(max_gen_cnt-1):

            # no need to go further if total max fitness is already 1
            total_max_fitness = self.find_max_fitness_and_child(i, True)
            if total_max_fitness['fitness'] >= 1:
                print "Success: "
                print total_max_fitness['child']
                break

            print "total max fitness"
            print total_max_fitness

            # survive half of the generation and pass them to the next gen
            self.survive_n_fittest(i, self.population_size / 2)


            # cross over the next generation half

            new_start_ind = self.cross_over_generation(i + 1, self.population_size / 2)
            # cross over twice to keep population constant
            self.cross_over_generation(i + 1, new_start_ind)

            # mutate the new ones only, not the old ones
            self.mutate_generation(i + 1, self.population_size / 2, self.population_size)


            # find fittest unmeasured child
            max_fitness = self.find_max_fitness_and_child(i + 1, False)

            print "unmeasured max fitness"
            print max_fitness

            if max_fitness['fitness'] >= 1:
                print "Success: "
                print max_fitness['child']
                break
            else:
                # update values for fittest unmeasured child through the experiment
                # self.score_handler.assign_percentages_for_all_patients(fittest_child)
                self.update_fitness(max_fitness['child'])
                self.score_handler.update_measured(max_fitness['child'])

                self.experiment_cnt += 1



    def animate(self, gen_cnt):
        """
        Run the evolution for n generations
        :param n:
        :return:
        """

        fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})

        rects = [Rectangle(xy=np.random.rand(2) * 10, width=0.5, height=0.5, angle=90) for i in
                 range(self.population_size)]

        def draw_generation(frame, self, rects):

            # for i in range(frame):
            i = frame
            for j in range(self.population_size):

                ind = j

                child = self.population[i][j][0:4]
                ax.add_artist(rects[ind])

                rects[i].set_clip_box(ax.bbox)
                alpha = self.population[i][j][3]/255

                rgb = self.population[i][j][0:3]/255
                rects[ind].set_alpha(alpha)
                rects[ind].set_facecolor(rgb)
                rects[ind].set_width(self.get_fitness_value(child))
                rects[ind].set_height(self.get_fitness_value(child))



            return ax

        ani = anim.FuncAnimation(fig, draw_generation, frames= gen_cnt, fargs=(self, rects), interval=25)
        # anim.FuncAnimation(fig, draw_generation, frames=gen_cnt, fargs=(self, rects), interval=50)

        # anim.FuncAnimation(fig, draw_generation, frames=10, repeat=False)

        ax.set_xlim(0, 10)
        ax.set_ylim(0, 10)


        plt.show()


#
gs = GASolver(MAX_GENERATIONS, POPULATION_SIZE, ANTIBODY_CNT, NC_CNT, C_CNT, MARKERS, C_THRESHOLD, CELL_CNT, CANCER_MU, CANCER_STD_DEV,
              NON_CANCER_MU, NON_CANCER_STD_DEV)

gs.run_simulation(MAX_GENERATIONS)

# gs.animate(gs.experiment_cnt)
