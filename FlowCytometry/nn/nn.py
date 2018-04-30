
# Importing Keras Sequential Model
from keras.models import Sequential
from keras.layers import Dense
import matplotlib.pyplot as plt
import numpy as np
from ga.scoreHandler import ScoreHandler
from util.staticMethods import StaticMethods
import operator

ANTIBODY_CNT = 64  # 16 # 64  # 240

NC_CNT = 5  # non-cancer patients
C_CNT = 5  # cancer patients


CELL_CNT = 100

DATA_SIZE = ANTIBODY_CNT / 4  # 40  # make this divisible by 4
MARKERS_LIST = [[0, 5, 10, 15]]
MU_LIST = [0.99]
STD_DEV_LIST = [0.01]

MUTATION_PROBABILITY = 0.1  # 0.03
CROSS_OVER_2_RATIO = 0.5

INPUT_DIM = 4


class NNSolver:

    def __init__(self, data_size, antibody_cnt, nc_cnt, c_cnt, markers_list, cell_cnt,
                 mu_list, sigma_list):

        self.data_size = data_size
        self.antibody_cnt = antibody_cnt
        self.nc_cnt = nc_cnt
        self.c_cnt = c_cnt
        self.data = np.full((self.data_size, INPUT_DIM), 0, dtype=np.int)
        self.output_list = np.full((self.data_size, 1), 0, dtype=np.float)

        self.score_handler = ScoreHandler(antibody_cnt, nc_cnt, c_cnt, markers_list, cell_cnt,
                                          mu_list, sigma_list)

        self.init_data()
        self.run_simulation()

        return

    def init_data(self):
        i = 0
        while i < ANTIBODY_CNT:
            ab_arr = np.array([i, i + 1, i + 2, i + 3])
            ab_arr = np.sort(ab_arr)
            self.data[i / INPUT_DIM] = ab_arr
            self.score_handler.update_measured(ab_arr)
            self.output_list[i / INPUT_DIM] = self.score_handler.compute_is_precise_for_ab_combination(ab_arr)
            i += 4

    def run_simulation(self):
        # Loading the data set (PIMA Diabetes Dataset)
        # dataset = np.loadtxt('datasets/pima-indians-diabetes.csv', delimiter=",")
        #
        # # Loading the input values to X and Label values Y using slicing.
        # X = dataset[:, 0:4]
        # Y = dataset[:, 4]

        # X = self.data[:, 0:INPUT_DIM]
        # Y = self.labels[:, 0:INPUT_DIM]
        #

        X = self.data
        Y = self.output_list

        print X
        print Y

        # Initializing the Sequential model from KERAS.
        model = Sequential()

        # Creating a 16 neuron hidden layer with Linear Rectified activation function.
        model.add(Dense(16, input_dim=4, kernel_initializer='uniform', activation='relu'))

        # Creating a 8 neuron hidden layer.
        model.add(Dense(8,  kernel_initializer='uniform', activation='relu'))

        # Adding a output layer.
        model.add(Dense(1,  kernel_initializer='uniform', activation='sigmoid'))

        # Compiling the model
        model.compile(loss='binary_crossentropy',
                      optimizer='adam', metrics=['accuracy'])
        # Fitting the modelz
        model.fit(X, Y, nb_epoch=150, batch_size=10)

        scores = model.evaluate(X, Y)

        print("%s: %.2f%%" % (model.metrics_names[1], scores[1] * 100))



nn = NNSolver( DATA_SIZE, ANTIBODY_CNT, NC_CNT, C_CNT, MARKERS_LIST, CELL_CNT, MU_LIST, STD_DEV_LIST)
