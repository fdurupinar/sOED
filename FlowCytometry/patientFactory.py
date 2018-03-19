from patient import Patient



class PatientFactory:

    def __init__(self, c_or_nc, ab_cnt,  markers, cnt, cell_cnt, c_mu, c_sigma, nc_mu, nc_sigma):
        """
        Create cnt patients of type <type>
        :param type: "c" or "nc"
        :param cnt: patient count to create
        """
        self.cell_cnt = cell_cnt
        self.ab_cnt = ab_cnt
        if c_or_nc == "c":
            self.mu = c_mu
            self.sigma = c_sigma
        else:
            self.mu = nc_mu
            self.sigma = nc_sigma

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