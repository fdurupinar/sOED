
# Importing Keras Sequential Model
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import LSTM
import numpy as np
from ga.scoreHandlerHypothetical import ScoreHandlerHypothetical
from itertools import combinations
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

from util.staticMethods import StaticMethods
import operator

ANTIBODY_CNT = 256  # 16 # 64  # 240

NC_CNT = 20  # non-cancer patients
C_CNT = 20  # cancer patients


CELL_CNT = 100

DATA_SIZE = ANTIBODY_CNT / 4  # 40  # make this divisible by 4
MARKERS_LIST = [[0, 5, 10, 15]]
MU_LIST = [0.9]
STD_DEV_LIST = [0.1]


INPUT_DIM = 4

NUM_EPOCHS = 150

class NNSolver:

    def __init__(self, data_size, antibody_cnt, nc_cnt, c_cnt, markers_list, cell_cnt,
                 mu_list, sigma_list):

        self.data_size = data_size
        self.antibody_cnt = antibody_cnt
        self.nc_cnt = nc_cnt
        self.c_cnt = c_cnt
        self.data = np.full((self.data_size, INPUT_DIM), 0, dtype=np.int)
        self.random_data = np.full((self.data_size, INPUT_DIM), 0, dtype=np.int)
        self.output_list = np.full((self.data_size, 1), 0, dtype=np.float)

        self.score_handler = ScoreHandlerHypothetical(antibody_cnt, nc_cnt, c_cnt, markers_list, cell_cnt,
                                                      mu_list, sigma_list)

        self.init_training_data()
        self.init_test_data()
        self.run_simulation()

        return

    def init_training_data(self):
        i = 0
        while i < ANTIBODY_CNT:
            ab_arr = np.array([i, i + 1, i + 2, i + 3])
            ab_arr = np.sort(ab_arr)
            self.data[i / INPUT_DIM] = ab_arr
            self.score_handler.update_measured(ab_arr)
            # self.output_list[i / INPUT_DIM] = self.score_handler.compute_max_precision_for_ab_combination(ab_arr)
            self.output_list[i / INPUT_DIM] = self.score_handler.compute_threshold_based_precision_for_ab_combination(ab_arr)
            # self.output_list[i / INPUT_DIM] = self.score_handler.compute_is_precise_for_ab_combination(ab_arr)
            i += 4



    def init_test_data(self):
        ab_names = np.arange(0, self.antibody_cnt)
        self.ab_combs = np.array([np.array(list(c)) for c in combinations(ab_names, 4)])  # combinations of 4
        # self.ab_combs = [list(c) for c in combinations(ab_names, 4)]  # combinations of 4

        # for comb in self.ab_combs:
        #     prec = self.score_handler.compute_threshold_based_precision_for_ab_combination(comb)
        #     print prec

    #     just to see the output

    def random_data(self):
        i = 0
        while i < ANTIBODY_CNT / 4:
            ab_arr = np.random.permutation(ANTIBODY_CNT)
            ab_arr = np.sort(ab_arr)

            while ab_arr in self.random_data:
                ab_arr = np.random.permutation(ANTIBODY_CNT)
                ab_arr = np.sort(ab_arr)

            self.random_data[i] = ab_arr
            i += 1

    def train_data(self, X, Y):
        # Initializing the Sequential model from KERAS.
        self.model = Sequential()

        # Creating a 16 neuron hidden layer with Linear Rectified activation function.
        self.model.add(Dense(16, input_dim=INPUT_DIM, kernel_initializer='uniform', activation='relu'))

        # Creating a 8 neuron hidden layer.
        self.model.add(Dense(8, kernel_initializer='uniform', activation='relu'))

        # Adding an output layer.
        self.model.add(Dense(1, kernel_initializer='uniform', activation='softmax'))

        # Compiling the model
        # This one is used for classification problems, i.e. not suitable for our precision-based problem
        # self.model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
        self.model.compile(loss='mean_squared_error', optimizer='adam', metrics=['accuracy'])

        # self.model.compile(loss='mean_squared_error', optimizer='adam', metrics=['accuracy'])
        # Fitting the modelz
        self.model.fit(X, Y, epochs=NUM_EPOCHS, batch_size=ANTIBODY_CNT/4)

        scores = self.model.evaluate(X, Y)

        print("\n%s: %.2f%%" % (self.model.metrics_names[1], scores[1] * 100))

    def run_simulation(self):

        # training dataGeneration
        X = self.data
        Y = self.output_list


        # print X
        print Y


        self.train_data(X, Y)

        # # calculate predictions
        # predictions = self.model.predict(self.ab_combs)
        # # round predictions
        # rounded = [round(x[0]) for x in predictions]
        # print(rounded)



        #

        # Generator





        # print("%s: %.2f%%" % (model_disc.metrics_names[1], scores[1] * 100))





nn = NNSolver( DATA_SIZE, ANTIBODY_CNT, NC_CNT, C_CNT, MARKERS_LIST, CELL_CNT, MU_LIST, STD_DEV_LIST)
