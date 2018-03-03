
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
import matplotlib.animation as anim
from scoreGenerator import ScoreGenerator

from itertools import combinations

ANTIBODY_CNT = 20 #242
MUTATION_PROBABILITY = 0.03
POPULATION_SIZE = 20  # make this divisible by 2
MAX_GENERATIONS = 100
NC_CNT = 20  # non-cancer patients
C_CNT = 20  # cancer patients

C_THRESHOLD = 0.4

class GASolver:

    def __init__(self):

        self.experiment_cnt = 0

        self.ab_combination_scores = np.zeros([ANTIBODY_CNT, ANTIBODY_CNT, ANTIBODY_CNT, ANTIBODY_CNT, 16])

        self.population = np.zeros([MAX_GENERATIONS, POPULATION_SIZE, 4])
        self.fitness = np.zeros([ANTIBODY_CNT, ANTIBODY_CNT, ANTIBODY_CNT, ANTIBODY_CNT])  # fitness of the combination

        self.current_generation = 0

        self.total_population = POPULATION_SIZE

        self.SG = ScoreGenerator(ANTIBODY_CNT, NC_CNT, C_CNT)

        self.create_random_population()

    def create_random_population(self):
        """
        Initial creation of a population of size POPULATION_SIZE
        :return:
        """

        # draw 4 random values from 242 antibodies n times
        # self.population = np.random.random_integers(ANTIBODY_CNT-1, size=(n, 4))

        for i in range(POPULATION_SIZE):
            ab_arr = np.random.permutation(np.arange(1, ANTIBODY_CNT))[0:4]

            # sortc
            ab_arr = np.sort(ab_arr)
            while ab_arr.tolist() in self.population[0].tolist():
                ab_arr = np.random.permutation(np.arange(1, ANTIBODY_CNT))[0:4]
                ab_arr = np.sort(ab_arr)

            self.population[0][i] = ab_arr

            # initial assignment of percentages for all the patients for the ab_arr and its combinations
            self.SG.assign_percentages_for_all_patients(ab_arr)


        # once the percentages are assigned, compute the precision scores deriving unknown values as multiplied
        # probabilities of known values

        for i in range(POPULATION_SIZE):
            self.update_fitness(self.population[0][i])

    def cross_over_generation(self, generation):
        """
        Cross over gene groups with the given ratio
        :param generation: generation number to reproduce
        :return:
        """

        # randomly divide these into two
        co_inds1 = np.random.permutation(POPULATION_SIZE)[0:POPULATION_SIZE/2]
        co_inds2 = np.setdiff1d(np.arange(POPULATION_SIZE), co_inds1)

        # cross over the next generation half
        last_ind = POPULATION_SIZE/2

        for i in range(len(co_inds1)):

            child = self.cross_over(self.population[generation][co_inds1[i]], self.population[generation][co_inds2[i]])

            if self.is_child_unique(generation + 1, child):
                self.population[generation + 1][last_ind] = child
                self.update_fitness(child)
            else:
                random_child = np.random.permutation(ANTIBODY_CNT-1)[0:4]
                while random_child.tolist() in self.population[generation+1].tolist():
                    random_child = np.random.permutation(ANTIBODY_CNT - 1)[0:4]
                    random_child = np.sort(random_child)
                self.population[generation + 1][last_ind] = random_child
                self.update_fitness(random_child)

            last_ind += 1
            self.total_population += 1  # increment total population




    def is_child_unique(self, generation, child):
        res = True
        for j in range(len(self.population[generation])):
            if child.tolist() in self.population[generation][j].tolist():
                res = False
        return res


    def cross_over(self, group1, group2):
        """
        Perform 2-point cross-over between two groups of antibodies of size 4
        :param group1:
        :param group2:
        :return:
        """

        child = np.zeros(4)

        inds1 = np.random.permutation(4)[0:2]  #select two indices
        inds2 = np.random.permutation(4)[0:2]  # select two indices

        ind = 0
        max_cnt = 2 # do this permutation max
        while (group1[inds1[0]] in group2.tolist() or group1[inds1[1]] in group2.tolist()) and ind < max_cnt:
            inds1 = np.random.permutation(4)[0:2]  # select two indices
            inds2 = np.random.permutation(4)[0:2]  # select two indices
            ind +=1

        if ind < max_cnt:
            # first half
            child[0] = group1[inds1[0]]
            child[1] = group1[inds1[1]]
            # second half
            child[2] = group2[inds2[0]]
            child[3] = group2[inds2[1]]
        else:
            child = np.random.permutation(ANTIBODY_CNT - 1)[0:4]

        child = np.sort(child)

        return child

    def mutate(self, child):
        """
        Randomly mutate one antibody in a group of 4
        :param child:
        :return:
        """
        mut_ind = np.random.random_integers(3)

        ab = np.random.random_integers(ANTIBODY_CNT-1)

        # keep on random selection if the current antibody already exists
        while ab in child:
            ab = np.random.random_integers(ANTIBODY_CNT -1)

        child[mut_ind] = ab
        child = np.sort(child)

        return child

    def increment_current_generation(self):
        """
        Increment current generation
        :return:
        """
        self.current_generation += 1

    def mutate_generation(self, generation):
        """Mutate some portion of the new generation"""

        for i in range(len(self.population[generation])/2, len(self.population[generation])):
            child = self.population[generation][i]

            mutation_chance = np.random.rand(1)[0]
            if mutation_chance < MUTATION_PROBABILITY:
                mutated_child = self.mutate(child)
                # if child exists, generate another child
                while self.is_child_unique(generation, mutated_child) == False:
                    mutated_child = self.mutate(mutated_child)
                self.population[generation][i] = mutated_child

    def update_fitness(self, child):
        """
        Compute the fitness of antibody quadruple based on existing values
        :param child:
        :return:
        """

        score = self.SG.compute_max_precision_for_ab_combination(child, C_THRESHOLD)


        self.fitness[child[0]][child[1]][child[2]][child[3]] = score

    def survive_n_fittest(self, generation, n):
        """Find and survive the n fittest combinations in the population, delete the rest"""

        fitness = np.zeros(len(self.population[generation]))
        for i in range(len(self.population[generation])):
            child = self.population[generation][i]
            fitness[i] = self.fitness[child[0]][child[1]][child[2]][child[3]]  # self.compute_fitness(child)

        fitness = np.sort(fitness)
        threshold = fitness[n-1]

        j = 0
        for i in range(len(self.population[generation])):
            child = self.population[generation][i]
            if self.fitness[child[0]][child[1]][child[2]][child[3]] > threshold:
                self.population[generation + 1][j] = child
                j += 1

        return j


    def run_simulation(self, max_gen_cnt):
        """
        Run evolution for max_gen_cnt iterations
        :param max_gen_cnt:
        :return:
        """




        for i in range(max_gen_cnt-1):

            prev_max_fitness = 1000

            iter_cnt = 0
            while True:
                max_fitness = -1
                fittest_child = self.population[i][0]
                for j in range(len(self.population[i])):
                    child = self.population[i][j]
                    fitness = self.fitness[child[0]][child[1]][child[2]][child[3]] # self.compute_fitness(self.population[i][j])

                    if fitness > max_fitness and fitness < prev_max_fitness:
                        max_fitness = fitness
                        fittest_child = child

                prev_max_fitness = max_fitness
                iter_cnt += 1
                # break the while loop if we find an unmeasured sequence
                if self.SG.measured[fittest_child[0]][fittest_child[1]][fittest_child[2]][fittest_child[3]] == 0 or \
                                iter_cnt >= len(self.population) :  # all of them have been measured
                    break

            prev_max_fitness = max_fitness

            print max_fitness
            print fittest_child

            # update values for fittest unmeasured child through the experiment
            self.SG.assign_percentages_for_all_patients(fittest_child)
            self.update_fitness(fittest_child)


            self.experiment_cnt += 1


            if max_fitness >= 1:
                break



            # survive half of the generation and pass them to the next gen
            self.survive_n_fittest(i, POPULATION_SIZE / 2)

            # cross over twice to keep population constant
            self.cross_over_generation(i)
            self.cross_over_generation(i)

            self.increment_current_generation()

            self.mutate_generation(i + 1)


    def animate(self, gen_cnt):
        """
        Run the evolution for n generations
        :param n:
        :return:
        """

        fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})


        rects = [Rectangle(xy=np.random.rand(2) * 10, width=0.5, height=0.5, angle=90) for i in
                 range(POPULATION_SIZE)]

        def draw_generation(frame, self, rects):

            # for i in range(frame):
            i = frame
            for j in range(POPULATION_SIZE):

                ind = j

                child = self.population[i][j][0:4]
                ax.add_artist(rects[ind])

                rects[i].set_clip_box(ax.bbox)
                alpha = self.population[i][j][3]/255

                rgb = self.population[i][j][0:3]/255
                rects[ind].set_alpha(alpha)
                rects[ind].set_facecolor(rgb)
                rects[ind].set_width(self.fitness[child[0]][child[1]][child[2]][child[3]])
                rects[ind].set_height(self.fitness[child[0]][child[1]][child[2]][child[3]])



            return ax

        ani = anim.FuncAnimation(fig, draw_generation, frames= gen_cnt, fargs=(self, rects), interval=25)
        # anim.FuncAnimation(fig, draw_generation, frames=gen_cnt, fargs=(self, rects), interval=50)

        # anim.FuncAnimation(fig, draw_generation, frames=10, repeat=False)

        ax.set_xlim(0, 10)
        ax.set_ylim(0, 10)


        plt.show()




gs = GASolver()
#
gs.run_simulation(MAX_GENERATIONS)
# gs.animate(gs.experiment_cnt)

# print(gs.population[0])

# gs.cross_over_population(0)

# gs.animate(25)
# gs.visualize_population()


# print(gs.population[1])
#
# #
# child = gs.cross_over([1,2,3,4], [5,6,7,8])
# print child
# child2 = gs.mutate(child)
# print child2

#
# print gs.population[0]
# gs.mutate_generation(0)
# print gs.population[0]