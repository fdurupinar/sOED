from patient import Patient
FILE_PATH = 'pdac.csv'

class DataSorter:

    def __init__(self, file_path):
        self.patients_c = []
        self.patients_nc = []

        self.read_data(file_path)

    def read_data(self, file_path):
        marker_name_list = []
        marker_id_list = []  # todo: will consider 4-plexes
        marker_count_benign = []
        marker_count_pdac = []
        marker_count_control = []


        with open(file_path, 'r') as fp:
            next(fp)  # skip first line

            ind = 0
            for line in fp:
                cols = line.split(',')
                marker_name_list.append(cols[0])
                marker_count_control.append(cols[3])
                marker_count_benign.append(cols[4])
                marker_count_pdac.append(cols[5])
                marker_id_list.append(ind)
                ind += 1


        p1_nc = Patient(False, marker_id_list, marker_count_benign)
        self.patients_nc.append(p1_nc)

        p1_c = Patient(True, marker_id_list, marker_count_pdac)
        self.patients_c.append(p1_c)
