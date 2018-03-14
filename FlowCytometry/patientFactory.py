from patient import Patient


CELL_CNT = 100

CANCER_MU = 0.8
CANCER_STD_DEV = 0.1

NON_CANCER_MU = 0.2
NON_CANCER_STD_DEV = 0.1


class PatientFactory:

    def __init__(self, c_or_nc, ab_cnt,  markers, cnt):
        """
        Create cnt patients of type <type>
        :param type: "c" or "nc"
        :param cnt: patient count to create
        """
        self.cell_cnt = CELL_CNT
        self.ab_cnt = ab_cnt
        if c_or_nc == "c":
            self.mu = CANCER_MU
            self.sigma = CANCER_STD_DEV
        else:
            self.mu = NON_CANCER_MU
            self.sigma = NON_CANCER_STD_DEV

        self.patients = self.create_patients(markers, cnt)

    def create_patients(self, markers, cnt):
        """
        Create cnt patients of type
        :param type:
        :param markers:
        :param cnt:
        :return:
        """

        patients = []

        for i in range(cnt):
            p = Patient(self.cell_cnt, self.ab_cnt, markers, self.mu, self.sigma)
            patients.append(p)

        return patients



# pfc = PatientFactory("c", [0, 1], 10)
# pfnc = PatientFactory("nc", [0, 1], 10)
#
# for p in pfc.patients:
#     for c in p.cells:
#         print c
#     print "c"
#
# for p in pfnc.patients:
#     for c in p.cells:
#         print c
#     print "nc"