from FlowCytometryTools import FCMeasurement
import os
import numpy as np
from pylab import *
from FlowCytometryTools import FCPlate
# Though the following import is not directly being used, it is required
# for 3D projection to work
from mpl_toolkits.mplot3d import Axes3D

from sklearn.cluster import KMeans

# import cytoflow as flow

path = '../data/'


# antibodies = {'fitc':0, 'apc':0, 'pe':0}
samples = np.zeros((60, 96))  # there will be max 60 patients


sample_dirs = os.listdir(path)

from FlowCytometryTools import ThresholdGate, PolyGate

#
sample_cnt = 0
for dir in sample_dirs:

    # get sample index
    ind_start = dir.find('#')
    if ind_start > -1:
        sample_ind = int(dir[ind_start+1:])
        sample_cnt += 1

        # read each antibody combination
        files = os.listdir(path+ '/' + dir)


        for datafile in files:
            # get datafile_id
            ab_ind = int(datafile[len(datafile) - 7: len(datafile) - 4])

            sample = FCMeasurement(ID=datafile, datafile=path+ '/' + dir + '/' + datafile)
            # sample = flow.Tube(file=path+ '/' + dir + '/' + datafile, conditions = {"Dox" : "float"})
            # samples[sample_ind] = sample


            # import_op = flow.ImportOp(conditions={"Dox": "float"},
            #                           tubes=[sample],
            #                           channels={'FITC-A': 'FITC-A', 'APC-A': 'APC-A', 'PE-A':'PE-A'})
            #
            #
            # k = flow.KMeansOp(name="KMeans",
            #                   channels=["FITC-A", "APC-A", "PE-A"],
            #                   scale={"FITC-A": "logicle","APC-A": "logicle", "PE-A": "logicle"},
            #                   num_clusters=8,
            #                   by=['Dox'])
            #
            #
            # ex = import_op.apply()
            #
            # k.estimate(ex)
            # ex2 = k.apply(ex)
            # k.default_view(yfacet="Dox").plot(ex2)
            #


            cells = sample.data[['FITC-A', 'APC-A', 'PE-A', 'SSC-A']].values

            fitc_gate = ThresholdGate(1000.0, ['FITC-A'], region='above')
            apc_gate = ThresholdGate(1000.0, ['APC-A'], region='above')

            figure();
            sample.plot(['FITC-A'], gates=[fitc_gate], bins=100);

            # fig = plt.figure()
            # ax = fig.add_subplot(111, projection='3d')

            # #
            # ax.scatter(cells[:,0], cells[:,1], cells[:,2])
            # #

            # plt.show()

            # tsample = sample.transform('hlog', channels=['FITC-A', 'APC-A', 'PE-A'], b=500.0)
            # figure()
            # tsample.plot(['FITC-A', 'APC-A'], kind = 'scatter', color='green', alpha=0.7, s=1);

            # print(len(sample))
            # print(len(cells))
            # print(cells[0])
            # print(cells[1])
            #
            # len(sample.data[''])
            # cell_cnt =
            #
            # samples[sample_ind][ab_ind] = np.full(cell_cnt, {'FITC':0, 'APC': 0, 'PE':0, 'SSC':0})

            # sample = FCMeasurement(ID='Test Sample', datafile=datafile)



        # samples[sample_ind] = np.zeros()
#
#
# # format is as follows:
# # sample[i] = [0: {'fitc':x, 'apc':y, 'pe':z}, 1: {'fitc':x, 'apc':y, 'pe':z}....]
#
#
# # datafiles = '../data/Sample #1/Specimen_001_F8_F08_068.fcs'
# #
# # sample = FCMeasurement(ID='Test Sample', datafile=datafile)
# #
# #
# # print(sample.channel_names)
# #
# # print(sample.data[['FITC-A', 'APC-A', 'PE-A']])
#
# # tsample = sample.transform('hlog', channels=['SSC-A', 'FITC-A', 'PE-A'], b=500.0)
# tsample = sample.transform('hlog', channels=['FITC-A', 'APC-A', 'PE-A'], b=500.0)