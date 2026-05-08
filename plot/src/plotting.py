import matplotlib.pyplot as plt

def hist(data, file_path, nbins=None, data_range=None, density=False):
    if (nbins is None and data_range is not None) or (nbins is not None and data_range is None):
        assert False
    
    plt.figure()
    if nbins is None:
        plt.hist(data, density=density)
    else:
        plt.hist(data, bins = [data_range[0] + (data_range[1] - data_range[0]) / (nbins + 1) * i for i in range(nbins + 1)], density=density)
    plt.savefig(file_path)
