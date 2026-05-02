import matplotlib.pyplot as plt

def extract_csv_cols(file_path, ncols):
    data = [[] for _ in range(ncols)]
    
    with open(file_path, mode='r') as file:
        for line in file.readlines():
            cols = line.split(", ")
            for i in range(ncols):
                data[i].append(float(cols[i]))

    return data

def plot_as_hist(data, file_path, nbins=None, data_range=None):
    if (nbins is None and data_range is not None) or (nbins is not None and data_range is None):
        assert False
    
    plt.figure()
    if nbins is None:
        plt.hist(data)
    else:
        plt.hist(data, bins = [data_range[0] + (data_range[1] - data_range[0]) / (nbins + 1) * i for i in range(nbins + 1)])
    plt.savefig(file_path)

if __name__ == '__main__':
    (m_values, cos_th_values, pT_values) = extract_csv_cols(
        '/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/data/events/events.lhe', 3)

    plot_as_hist(
        m_values, 
        '/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/plot/m_hist.pdf',
        nbins=20,
        data_range=(20, 200)
    )

    plot_as_hist(
        cos_th_values, 
        '/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/plot/cos_th_hist.pdf',
        nbins=20,
        data_range=(-1, 1)
    )

    plot_as_hist(
        pT_values, 
        '/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/plot/first_em_pT_hist.pdf',
        nbins=20,
        data_range=(0, 150)
    )
