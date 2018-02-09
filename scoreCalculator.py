from scipy.stats import norm
import matplotlib.pyplot as plt
import numpy as np


ANTIBODY_CNT = 20

#A,B,C,D combinations [-,-,-,-],[-,-,-,+], [-,-,+,-]. [-,-,+,+], [-,+,-,-]. [-,+,-,+]. [-,+,+,-], [-,+,+,+].
#[+,-,-,-],[+,-,-,+], [+,-,+,-]. [+,-,+,+], [+,+,-,-]. [+,+,-,+]. [+,+,+,-], [+,+,+,+]


class ScoreCalculator:
    def __init__(self):
        # fig, ax = plt.subplots(1, 1)

        self.mu, self.sigma = 0, 0.5  # mean and standard deviation
        # s = np.random.normal(self.mu, self.sigma, 1000) # draw 1000 samples
        #
        # count, bins, ignored = plt.hist(s, 50, normed=True)
        #
        # plt.plot(bins, 1 / (self.sigma * np.sqrt(2 * np.pi)) * np.exp(- (bins - self.mu) ** 2 / (2 * self.sigma ** 2)),
        #          linewidth = 2, color = 'r')
        # plt.show()
        self.antibodies = np.arange(1, ANTIBODY_CNT)
        self.ab_combination_scores = np.zeros([ANTIBODY_CNT, ANTIBODY_CNT, ANTIBODY_CNT, ANTIBODY_CNT, 16])
        self.ab_precision_values = np.zeros([ANTIBODY_CNT, ANTIBODY_CNT, ANTIBODY_CNT, ANTIBODY_CNT, 16])


    # these are randomly drawn from a gaussian distribution for now, but will be obtained from the experiments
    def drawPrecisionScores(self, f_known):
        """Draw 16 values from normal distribution excluding known values"""

        f = np.random.normal(self.mu, self.sigma, 16)

        # TODO: do more efficiently
        # replace the already known values with f_known
        for i in range(0, 16):
            if f_known[i] != 0:
                f[i] = f_known[i]

        return f

    def calculateScore(self, f, c_known):
        """calculate the coefficients for 16 precision values that gives 99% recognition"""

        # compute the precision value based on existing values
        required_prec = 0.99 - np.dot(f, c_known)

        f_inv = 1/f/16

        c = f_inv * required_prec

        # TODO: do more efficiently
        # replace already known values with c_known
        for i in range(0,16):
            if c_known[i] != 0:
                c[i] = c_known[i]

        return c


    def fillDynamicTable(self, exp_cnt):
        """exp_cnt = experiment count"""

        for i in range(0, exp_cnt):
            # draw 4 random values from 242 antibodies and fill in the table related to their combinations
            ab = np.zeros(4)
            for j in range(0,4):
                val = np.random.random_integers(ANTIBODY_CNT-1-j)
                while val in ab:
                    val = np.random.random_integers(ANTIBODY_CNT-1-j)
                ab[j] = val

            # sort ab so that all the elements are in order for the antibody combinations
            ab = np.sort(ab)
            # print(str(ab[0]) + " " + str(ab[1]) + " " + str(ab[2]) + " " + str(ab[3]))

            c_known = self.ab_combination_scores[ab[0]][ab[1]][ab[2]][ab[3]]

            f_known = self.ab_precision_values[ab[0]][ab[1]][ab[2]][ab[3]]

            # generate scores
            f = self.drawPrecisionScores(f_known)
            c = self.calculateScore(f, c_known)

            self.ab_combination_scores[ab[0]][ab[1]][ab[2]][ab[3]] = c
            self.ab_precision_values[ab[0]][ab[1]][ab[2]][ab[3]] = f

            # print(self.ab_combination_scores[ab[0]][ab[1]][ab[2]][ab[3]])




sc = ScoreCalculator()

sc.fillDynamicTable(260)


# for i in range(0, ANTIBODY_CNT):
#     for j in range(0, ANTIBODY_CNT):
#         for k in range(0, ANTIBODY_CNT):
#             for l in range(0, ANTIBODY_CNT):
#                 # if np.linalg.norm(sc.ab_combination_scores[i][j][k][l]) != 0:
#                     print(sc.ab_combination_scores[i][j][k][l])

