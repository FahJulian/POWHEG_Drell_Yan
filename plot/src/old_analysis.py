import plotting

def extract_csv_cols(file_path, ncols):
    data = [[] for _ in range(ncols)]
    
    with open(file_path, mode='r') as file:
        for line in file.readlines():
            cols = line.split(", ")
            for i in range(ncols):
                data[i].append(float(cols[i]))

    return data

def old_analysis():
    (m_values, cos_th_values, pT_values) = extract_csv_cols(
        '/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/data/events/events_100k.lhe', 3)

    plotting.hist(
        m_values, 
        '/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/plot/m_hist.pdf',
        nbins=20,
        data_range=(75, 105),
        density=True
    )

    plotting.hist(
        cos_th_values, 
        '/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/plot/cos_th_hist.pdf',
        nbins=20,
        data_range=(-1, 1),
        density=True
    )

    plotting.hist(
        pT_values, 
        '/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/plot/first_em_pT_hist.pdf',
        nbins=20,
        data_range=(0, 150),
        density=True
    )
