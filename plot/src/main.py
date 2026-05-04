import matplotlib.pyplot as plt
import numpy as np

def extract_csv_cols(file_path, ncols):
    data = [[] for _ in range(ncols)]
    
    with open(file_path, mode='r') as file:
        for line in file.readlines():
            cols = line.split(", ")
            for i in range(ncols):
                data[i].append(float(cols[i]))

    return data

def remove_whitespace(line):
    line = line.replace('\n', '')
    line = line.strip()

    while '  ' in line:
        line = line.replace('  ', ' ')

    return line

def extract_pT(lepton_momenta):
    if len(lepton_momenta) != 2:
        assert(False)

    total_momentum = lepton_momenta[0] + lepton_momenta[1]
    pT2 = total_momentum[1]**2 + total_momentum[2]**2

    assert(pT2 > -1e-12)
    return np.sqrt(max(pT2, 0.0))

def extract_from_lhe(file_path):
    in_event = False
    event_line = 1
    lepton_momenta = []
    boson_mass = None

    m_values = []
    pT_values = []
    cos_th_values = []

    with open(file_path, mode='r') as file:
        for line in file.readlines():
            line = remove_whitespace(line)

            if line.startswith('#'):
                continue

            if in_event:
                if event_line != 1:
                    if line.startswith('11') or line.startswith('-11'):
                        cols = line.split(' ')
                        e = float(cols[9])
                        pX = float(cols[6])
                        pY = float(cols[7])
                        pZ = float(cols[8])
                        lepton_momenta.append(np.array([e, pX, pY, pZ]))
                    elif line.startswith('23'):
                        cols = line.split(' ')
                        boson_mass = float(cols[10])

                event_line += 1

            if line == '<event>': 
                in_event = True

                lepton_momenta = []
                boson_mass = None
                event_line = 1
            elif line == '</event>':
                pT_values.append(extract_pT(lepton_momenta))
                m_values.append(boson_mass)
                in_event = False
            
            if len(pT_values) % 1000 == 0:
                print(str(len(pT_values)) + ' Events read')

    return m_values, cos_th_values, pT_values

def plot_as_hist(data, file_path, nbins=None, data_range=None, density=False):
    if (nbins is None and data_range is not None) or (nbins is not None and data_range is None):
        assert False
    
    plt.figure()
    if nbins is None:
        plt.hist(data, density=density)
    else:
        plt.hist(data, bins = [data_range[0] + (data_range[1] - data_range[0]) / (nbins + 1) * i for i in range(nbins + 1)], density=density)
    plt.savefig(file_path)

if __name__ == '__main__':
    m_values, cos_th_values, pT_values = extract_from_lhe(
        '/home/julian/documents/uni/master/master_thesis/learning/powheg/out/z/pwgevents.lhe'
    )

    plot_as_hist(
        m_values, 
        '/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/plot/m_hist_pwhg.pdf',
        nbins=20,
        data_range=(20, 200),
        density=True
    )

    plot_as_hist(
        cos_th_values, 
        '/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/plot/cos_th_hist_pwhg.pdf',
        nbins=20,
        data_range=(-1, 1),
        density=True
    )

    plot_as_hist(
        pT_values, 
        '/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/plot/first_em_pT_hist_pwhg.pdf',
        nbins=20,
        data_range=(0, 150),
        density=True
    )

    (m_values, cos_th_values, pT_values) = extract_csv_cols(
        '/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/data/events/events.lhe', 3)

    plot_as_hist(
        m_values, 
        '/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/plot/m_hist.pdf',
        nbins=20,
        data_range=(20, 200),
        density=True
    )

    plot_as_hist(
        cos_th_values, 
        '/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/plot/cos_th_hist.pdf',
        nbins=20,
        data_range=(-1, 1),
        density=True
    )

    plot_as_hist(
        pT_values, 
        '/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/plot/first_em_pT_hist.pdf',
        nbins=20,
        data_range=(0, 150),
        density=True
    )

