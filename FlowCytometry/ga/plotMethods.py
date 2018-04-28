import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
import matplotlib.animation as anim


class PlotMethods:

    def __init__(self, ga_solver):
        self.ga_solver = ga_solver

    def _ab_to_xy(self, ab_arr):
        """
        Generate a unique tuple to encode 4 antibodies
        :param ab_arr:
        :return: [x, y] representing antibodies' unique position in 2d space
        """

        x = ab_arr[0] + ab_arr[1] * self.ga_solver.antibody_cnt
        y = ab_arr[2] + ab_arr[3] * self.ga_solver.antibody_cnt

        return [x, y]

    def visualize_generation(self, generation):
        """
        Draw each ab combination color and size coded
        :param generation:
        :param ax:
        :return:
        """

        fig, ax = plt.subplots()

        # ax.clear()
        for ind in range(len(self.ga_solver.population[generation])):

            child = self.ga_solver.population[generation][ind]

            xy = self._ab_to_xy(child)
            rect = Rectangle(xy=xy, width=0.5, height=0.5, angle=90)
            ax.add_artist(rect)

            rect.set_clip_box(ax.bbox)
            # alpha = float(child[3]) / ANTIBODY_CNT
            # rgb = [float(child[0])/ANTIBODY_CNT, float(child[1])/ANTIBODY_CNT, float(child[2])/ANTIBODY_CNT/ ANTIBODY_CNT]
            #
            # rect.set_alpha(alpha)
            # rect.set_facecolor(rgb)

            fitness = self.ga_solver.get_fitness_value(child)

            if fitness > 0.5:
                rect.set_facecolor('r')
            else:
                rect.set_facecolor('b')

            rect.set_width(fitness * 100)
            rect.set_height(fitness * 100)

        ax.set_xlim(0, self.ga_solver.antibody_cnt * self.ga_solver.antibody_cnt + 10)
        ax.set_ylim(0, self.ga_solver.antibody_cnt * self.ga_solver.antibody_cnt + 10)

        plt.show()
        plt.pause(1)
        plt.close()

    def plot_generation(self, generation):
        """
        Color code each antibody combination as RGBA
        :param generation:
        :return:
        """
        fig, ax = plt.subplots(figsize=(10, 10))
        x_data = np.arange(0, self.ga_solver.population_size)

        x_labels = []
        y_data = []
        color_data = []
        for i in range(self.ga_solver.population_size):
            child = self.ga_solver.population[generation][i]
            fitness = self.ga_solver.get_fitness_value(child)

            y_data.append(fitness)
            x_labels.append(str(child[0]) + "-" + str(child[1]) + "-" + str(child[2]) + "-" + str(child[3]))
            color_data.append(
                [float(child[0]) / self.ga_solver.antibody_cnt , float(child[1]) / self.ga_solver.antibody_cnt , float(child[2]) / self.ga_solver.antibody_cnt,
                 float(child[3]) / self.ga_solver.antibody_cnt ])

        plt.bar(x_data, y_data, color=color_data, tick_label=x_labels)
        plt.xticks(rotation=90)
        plt.show()
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        plt.pause(1)
        plt.close()

    def animate(self, gen_cnt):
        """
        Run the evolution for n generations
        :param n:
        :return:
        """

        fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})

        rects = [Rectangle(xy=np.random.rand(2) * 10, width=0.5, height=0.5, angle=90) for i in
                 range(self.ga_solver.population_size)]

        def draw_generation(frame, self, rects):
            # for i in range(frame):
            i = frame
            for j in range(self.ga_solver.population_size):
                ind = j

                child = self.ga_solver.population[i][j][0:4]
                ax.add_artist(rects[ind])

                rects[i].set_clip_box(ax.bbox)
                alpha = self.ga_solver.population[i][j][3] / self.ga_solver.antibody_cnt

                rgb = self.ga_solver.population[i][j][0:3] / self.ga_solver.antibody_cnt
                rects[ind].set_alpha(alpha)
                rects[ind].set_facecolor(rgb)
                rects[ind].set_width(self.ga_solver.get_fitness_value(child))
                rects[ind].set_height(self.ga_solver.get_fitness_value(child))

            return ax

        ani = anim.FuncAnimation(fig, draw_generation, frames=gen_cnt, fargs=(self, rects), interval=25)

        ax.set_xlim(0, 10)
        ax.set_ylim(0, 10)

        plt.show()

    # def plot_input(self):
    #     """
    #     Plots the initial data with cells and ab distributions
    #     :return:
    #     """
    #     x_data = []
    #     y_data = []
    #
    #     # PLOT normal patients
    #     for nc in self.score_handler.patients_nc:
    #         for p in self.population[0]: # first generation
    #             child = p[0:4]
    #
    #
    #
    #     for i in range(self.ab_cnt):
    #
    #         x_data.append(i)
    #         y_data.append(self.get_marker_ratio([i]))
    #
    #     plt.plot(x_data, y_data)
    #
    #     plt.xlabel('Antibody index')
    #     plt.ylabel('Marker ratio')
    #     plt.grid(True)
    #
    #     plt.show()


#
