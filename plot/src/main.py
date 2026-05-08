import matplotlib.pyplot as plt
import numpy as np

import plotting
import lhe_reading

OUT_FOLDER = '/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/plot/'

def extract_pt_values(data : lhe_reading.Data):
    px_values = np.array([p[1] for p in data.momenta_boson])
    py_values = np.array([p[2] for p in data.momenta_boson])

    return np.sqrt(px_values * px_values + py_values * py_values)

if __name__ == '__main__':
    data = lhe_reading.read_lhe('/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/data/events/events.lhe')
    data_powheg = lhe_reading.read_lhe('/home/julian/documents/uni/master/master_thesis/learning/powheg/out/z/pwgevents.lhe')

    plotting.hist(data_powheg.masses_boson, OUT_FOLDER + 'm_hist_pwhg.pdf', nbins=20, data_range=(75, 105), density=True)
    plotting.hist(extract_pt_values(data_powheg), OUT_FOLDER + 'first_em_pT_hist_pwhg.pdf', nbins=20, data_range=(0, 150), density=True)
    # plotting.plot_as_hist(cos_th_values, OUT_FOLDER + 'cos_th_hist_pwhg.pdf', nbins=20, data_range=(-1, 1), density=True)

    plotting.hist(data.masses_boson, OUT_FOLDER + 'm_hist.pdf', nbins=20, data_range=(75, 105), density=True)
    plotting.hist(extract_pt_values(data), OUT_FOLDER + 'first_em_pT_hist.pdf', nbins=20, data_range=(0, 150), density=True)
    # plotting.plot_as_hist(cos_th_values, OUT_FOLDER + 'cos_th_hist_pwhg.pdf', nbins=20, data_range=(-1, 1), density=True)
