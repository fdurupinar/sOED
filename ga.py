
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
import matplotlib.animation as anim

ANTIBODY_CNT = 242
GENERATIONS = 1000
MUTATION_PROBABILITY = 0.03
POPULATION_SIZE = 20
CROSS_OVER_PERCENT = 0.5
MAX_GENERATIONS = 20
FITNESS_THRESHOLD = 0.8

class GASolver:

    def __init__(self):

        self.mu, self.sigma = 0, 0.5  # mean and standard deviation

        self.antibodies = np.arange(1, ANTIBODY_CNT)
        self.ab_combination_scores = np.zeros([ANTIBODY_CNT, ANTIBODY_CNT, ANTIBODY_CNT, ANTIBODY_CNT, 16])
        self.ab_precision_values = np.zeros([ANTIBODY_CNT, ANTIBODY_CNT, ANTIBODY_CNT, ANTIBODY_CNT, 16])
        self.population = np.zeros([MAX_GENERATIONS, POPULATION_SIZE, 4])
        self.current_generation = 0
        self.create_random_population()
        self.total_population = POPULATION_SIZE

        # these are randomly drawn from a gaussian distribution for now, but will be obtained from the experiments
    def draw_precision_scores(self, n):
        """
        Draw 16 values from normal distribution n times
        :param n: sample count
        :return:
        """

        # for each sample
        f = np.zeros([n, 16])
        for i in range(n):
            f[i] = np.random.normal(self.mu, self.sigma, 16)

        return f

    def create_random_population(self):
        """
        Initial creation of a population of size POPULATION_SIZE
        :return:
        """

        # draw 4 random values from 242 antibodies n times
        # self.population = np.random.random_integers(ANTIBODY_CNT-1, size=(n, 4))

        for i in range(POPULATION_SIZE):
            ab_arr = np.random.permutation(ANTIBODY_CNT-1)[0:4]

            # sort
            ab_arr = np.sort(ab_arr)
            while ab_arr in self.population[0]:
                ab_arr = np.random.permutation(ANTIBODY_CNT - 1)[0:4]
                ab_arr = np.sort(ab_arr)

            self.population[0][i] = ab_arr

    def cross_over_generation(self, generation):
        """
        Cross over gene groups with the given ratio
        :param generation: generation number to reproduce
        :return:
        """

        # divide to half
        co_cnt = int(len(self.population[generation]) * CROSS_OVER_PERCENT / 2)


        # indices of population members who will cross over
        co_inds1 = np.random.random_integers(len(self.population[generation])-1, size=co_cnt)
        co_inds2 = np.random.random_integers(len(self.population[generation])-1, size=co_cnt)

        for i in range(co_cnt):
            child = self.cross_over(self.population[generation][co_inds1[i]], self.population[generation][co_inds2[i]])

            if self.is_child_unique(child):
                np.append(self.population[generation + 1], [child])
                self.total_population += 1  # increment total population



    def is_child_unique(self, child):
        for i in range(len(self.population)):
            for j in range(len(self.population[i])):
                if child in self.population[i][j]:
                    return False
        return True


    def cross_over(self, group1, group2):
        """
        Perform 2-point cross-over between two groups of antibodies of size 4
        :param group1:
        :param group2:
        :return:
        """

        inds1 = np.random.permutation(4)[0:2]  #select two indices
        inds2 = np.random.permutation(4)[0:2]  # select two indices

        child = np.zeros(4)
        # first half
        child[0] = group1[inds1[0]]
        child[1] = group1[inds1[1]]
        # second half
        child[2] = group2[inds2[0]]
        child[3] = group2[inds2[1]]

        child = np.sort(child)

        return child

    def mutate(self, child):
        """
        Randomly mutate one antibody in a group of 4
        :param child:
        :return:
        """
        mut_ind = np.random.random_integers(3)

        ab = np.random.random_integers(ANTIBODY_CNT - 1)

        # keep on random selection if the current antibody already exists
        while ab in child:
            ab = np.random.random_integers(ANTIBODY_CNT - 1)

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

        for child in self.population[generation]:
            mutation_chance = np.random.rand(1)[0]
            if mutation_chance < MUTATION_PROBABILITY:
                child = self.mutate(child)

                self.delete_child(generation, child)


    def delete_child(self, generation, child):
        """ Delete child if it exists"""

        ind = np.where(self.population[generation] == child)
        if ind and ind[0] and ind[0][0]:
            np.delete(self.population[generation], ind[0][0])

    def eliminate_weaklings(self):
        """
        Delete the unfit ones
        :return:
        """
        for i in range(len(self.population)):
            for child in self.population[i]:
                fitness = self.compute_fitness(child)
                if fitness < FITNESS_THRESHOLD:
                    self.delete_child(i, child)

    def compute_fitness(self, child):
        """
        Compute the fitness of antibody quadruple
        :param child:
        :return:
        """


    def run_simulation(self, max_gen_cnt):
        """
        Run evolution for max_gen_cnt iterations
        :param max_gen_cnt:
        :return:
        """

        for i in range(max_gen_cnt):
            self.cross_over_generation(i)
            self.increment_current_generation()
            self.mutate_generation(i + 1)
            self.eliminate_weaklings()


    def animate(self, gen_cnt):
        """
        Run the evolution for n generations
        :param n:
        :return:
        """


        fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})

        #
        # for i in range(n):
        #     self.cross_over_population(i)

        rects = [Rectangle(xy=np.random.rand(2) * 10, width=0.5, height=0.5, angle=90) for i in
                 range(self.total_population)]

        def draw_generation(frame, self, rects):
            # line.set_ydata(np.sin(x + frame / 10.0))  # update the data
            # return line,

            rects = [Rectangle(xy=np.random.rand(2) * 10, width=0.5, height=0.5, angle=90) for i in
                     range(self.total_population)]

            print(frame)

            for i in range(0, frame):
                for j in range(0, len(self.population[i])):

                    ind = i * len(self.population[i]) + j

                    ax.add_artist(rects[ind])

                    rects[i].set_clip_box(ax.bbox)
                    alpha = self.population[i][j][3]/255

                    rgb = self.population[i][j][0:3]/255
                    rects[ind].set_alpha(alpha)
                    rects[ind].set_facecolor(rgb)


            self.cross_over_generation(frame)

            return ax

        ani = anim.FuncAnimation(fig, draw_generation, frames= gen_cnt, fargs=(self, rects), interval=25)
        # anim.FuncAnimation(fig, draw_generation, frames=gen_cnt, fargs=(self, rects), interval=50)

        # anim.FuncAnimation(fig, draw_generation, frames=10, repeat=False)

        ax.set_xlim(0, 10)
        ax.set_ylim(0, 10)


        plt.show()



            # self.visualize_population()





    def visualize_population(self):
        """
        Visualize each member as a rgba-coded rectangle
        :param current_generation: Generation number to visualize up to
        :return:
        """



        fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})


        rects = [Rectangle(xy=np.random.rand(2) * 10, width=0.5, height=0.5, angle=90) for i in
                 range(self.total_population)]

        for i in range(0, self.current_generation):
            for j in range(0, len(self.population[i])):
                ind = i * len(self.population[i]) + j
                ax.add_artist(rects[ind])
                rects[i].set_clip_box(ax.bbox)
                alpha = self.population[i][j][3]/255

                rgb = self.population[i][j][0:3]/255
                rects[ind].set_alpha(alpha)
                rects[ind].set_facecolor(rgb)

        # anim.FuncAnimation(fig, draw_generation, frames=10, repeat=False)

        ax.set_xlim(0, 10)
        ax.set_ylim(0, 10)


        plt.show()





gs = GASolver()
gs.run_simulation(10)

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