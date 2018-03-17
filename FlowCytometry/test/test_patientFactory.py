from unittest import TestCase
from patientFactory import PatientFactory


cell_cnt = 100

c_mu = 0.7
c_sigma = 0.1

nc_mu = 0.3
nc_sigma = 0.1

ab_cnt = 5
cell_cnt = 10
markers = [1, 2]

class TestPatientFactory(TestCase):

    def test_init(self):
        pfnc = PatientFactory("nc", ab_cnt, markers, 5, cell_cnt, c_mu, c_sigma, nc_mu, nc_sigma)
        pfc = PatientFactory("c", ab_cnt, markers, 5, cell_cnt, c_mu, c_sigma, nc_mu, nc_sigma)

        self.assertEqual(pfnc.mu, nc_mu)
        self.assertEqual(pfnc.sigma, nc_sigma)

        self.assertEqual(pfc.mu, c_mu)
        self.assertEqual(pfc.sigma, c_sigma)

        self.assertEqual(len(pfnc.patients), 5)
        self.assertEqual(len(pfc.patients), 5)

        self.assertEqual(len(pfnc.patients[0].cells), cell_cnt)
        self.assertEqual(len(pfc.patients[0].cells), cell_cnt)
