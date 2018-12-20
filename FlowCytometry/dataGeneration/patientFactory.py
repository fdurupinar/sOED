from dataGeneration.patient import Patient


class PatientFactory:

    def __init__(self, patient_type, cnt, ab_cnt, markers_list,  cell_cnt, mu_list, sigma_list):
        """
        Create cnt patients of type <type>
        :param type: "c" or "nc"
        :param cnt: patient count to create
        """

        self.patients = []

        for i in range(cnt):
            p = Patient(patient_type, cell_cnt, ab_cnt, markers_list, mu_list, sigma_list)
            self.patients.append(p)

