
# Importing Keras Sequential Model
from keras.models import Sequential
from keras.layers import Dense
import matplotlib.pyplot as plt
import numpy as np
from util.staticMethods import StaticMethods
import operator

ANTIBODY_CNT = 64  # 16 # 64  # 240

NC_CNT = 5  # non-cancer patients
C_CNT = 5  # cancer patients


CELL_CNT = 100


MARKERS_LIST = [[0, 5, 10, 15]]
MU_LIST = [0.99]
STD_DEV_LIST = [0.01]

MUTATION_PROBABILITY = 0.1  # 0.03
CROSS_OVER_2_RATIO = 0.5
VISUALIZE_POPULATION = False


class NNSolver:

    def __init__(self, antibody_cnt, nc_cnt, c_cnt, markers_list,  cell_cnt,
                 mu_list, sigma_list):


        return

    def run_simulation(self):
        # Loading the data set (PIMA Diabetes Dataset)
        dataset = np.loadtxt('datasets/pima-indians-diabetes.csv', delimiter=",")

        # Loading the input values to X and Label values Y using slicing.
        X = dataset[:, 0:4]
        Y = dataset[:, 4]

        # Initializing the Sequential model from KERAS.
        model = Sequential()

        # Creating a 16 neuron hidden layer with Linear Rectified activation function.
        model.add(Dense(16, input_dim=4, init='uniform', activation='relu'))

        # Creating a 8 neuron hidden layer.
        model.add(Dense(8, init='uniform', activation='relu'))

        # Adding a output layer.
        model.add(Dense(1, init='uniform', activation='sigmoid'))

        # Compiling the model
        model.compile(loss='binary_crossentropy',
                      optimizer='adam', metrics=['accuracy'])
        # Fitting the model
        model.fit(X, Y, nb_epoch=150, batch_size=10)

        scores = model.evaluate(X, Y)

        print("%s: %.2f%%" % (model.metrics_names[1], scores[1] * 100))



