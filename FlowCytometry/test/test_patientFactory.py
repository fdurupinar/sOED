from unittest import TestCase
from patientFactory import PatientFactory
from patientFactory import CANCER_MU
from patientFactory import CANCER_STD_DEV
from patientFactory import NON_CANCER_MU
from patientFactory import NON_CANCER_STD_DEV
from patientFactory import CELL_CNT

ab_cnt = 5
cell_cnt = 10
markers = [1, 2]

class TestPatientFactory(TestCase):

    def test_init(self):
        pfnc = PatientFactory("nc", ab_cnt, markers, 5)
        pfc = PatientFactory("c", ab_cnt, markers, 5)

        self.assertEqual(pfnc.mu, NON_CANCER_MU)
        self.assertEqual(pfnc.sigma, NON_CANCER_STD_DEV)

        self.assertEqual(pfc.mu, CANCER_MU)
        self.assertEqual(pfc.sigma, CANCER_STD_DEV)

        self.assertEqual(len(pfnc.patients), 5)
        self.assertEqual(len(pfc.patients), 5)

        self.assertEqual(len(pfnc.patients[0].cells), CELL_CNT)
        self.assertEqual(len(pfc.patients[0].cells), CELL_CNT)
