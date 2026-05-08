import matplotlib.pyplot as plt
import numpy as np

import plotting
import lhe_reading

OUT_FOLDER = '/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/plot/100k_events/'

def extract_pt_values(momenta : np.array):
    px_values = np.array([p[1] for p in momenta])
    py_values = np.array([p[2] for p in momenta])

    return np.sqrt(px_values * px_values + py_values * py_values)

if __name__ == '__main__':
    data = lhe_reading.read_lhe('/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/data/events/events_10k.lhe')
    data_powheg = lhe_reading.read_lhe('/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/data/events/pwgevents_10k.lhe')

    plotting.hist(data_powheg.masses_boson, OUT_FOLDER + 'm_hist_pwhg.pdf', nbins=20, data_range=(75, 105), density=True)
    plotting.hist(extract_pt_values(data_powheg.momenta_boson), OUT_FOLDER + 'boson_pt_hist_pwhg.pdf', nbins=20, data_range=(0, 150), density=True)
    plotting.hist(extract_pt_values(data_powheg.momenta_l_minus), OUT_FOLDER + 'lepton1_pt_hist_pwhg.pdf', nbins=20, data_range=(0, 55), density=True)
    plotting.hist(extract_pt_values(data_powheg.momenta_l_plus), OUT_FOLDER + 'lepton2_pt_hist_pwhg.pdf', nbins=20, data_range=(0, 55), density=True)
    # plotting.plot_as_hist(cos_th_values, OUT_FOLDER + 'cos_th_hist_pwhg.pdf', nbins=20, data_range=(-1, 1), density=True)

    plotting.hist(data.masses_boson, OUT_FOLDER + 'm_hist.pdf', nbins=20, data_range=(75, 105), density=True)
    plotting.hist(extract_pt_values(data.momenta_boson), OUT_FOLDER + 'boson_pt_hist.pdf', nbins=20, data_range=(0, 150), density=True)
    plotting.hist(extract_pt_values(data.momenta_l_minus), OUT_FOLDER + 'lepton1_pt_hist.pdf', nbins=20, data_range=(0, 55), density=True)
    plotting.hist(extract_pt_values(data.momenta_l_plus), OUT_FOLDER + 'lepton2_pt_hist.pdf', nbins=20, data_range=(0, 55), density=True)
    # plotting.plot_as_hist(cos_th_values, OUT_FOLDER + 'cos_th_hist_pwhg.pdf', nbins=20, data_range=(-1, 1), density=True)
