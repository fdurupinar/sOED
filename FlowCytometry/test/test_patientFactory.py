from unittest import TestCase
from patientFactory import PatientFactory


cell_cnt = 100

c_mu_list = [0.7, 0.5]
c_sigma_list = [0.1, 0.1]

nc_mu_list = [0.3, 0.2]
nc_sigma_list = [0.1, 0.1]

ab_cnt = 5
cell_cnt = 10
markers_list = [[1, 2], [3]]


class TestPatientFactory(TestCase):

    def test_init(self):
        pfnc = PatientFactory(5, ab_cnt, markers_list, cell_cnt, c_mu_list, c_sigma_list)
        pfc = PatientFactory(5, ab_cnt, markers_list,  cell_cnt, nc_mu_list, nc_sigma_list)

        self.assertEqual(len(pfnc.patients), 5)
        self.assertEqual(len(pfc.patients), 5)

        self.assertEqual(len(pfnc.patients[0].cells), cell_cnt)
        self.assertEqual(len(pfc.patients[0].cells), cell_cnt)
