import numpy as np
from tqdm import tqdm
import os
import h5py
import argparse

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit

def gauss(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

# Create a 2d histogram and save it to the folder `direc` with the name `filename`.
def occupancy(x_coords, y_coords, direc, filename):
    # Definition of the plot size
    fig_width, fig_height = 7, 6  # in inches

    # Relative font size
    font_size = fig_height * 2

    # Set font size
    plt.rcParams.update({
        'font.size': font_size,
        'axes.titlesize': font_size,
        'axes.labelsize': font_size,
        'xtick.labelsize': font_size,
        'ytick.labelsize': font_size,
        'legend.fontsize': font_size,
    })

    p0_x = [1000, np.mean(x_coords), 10]
    hist_x, bin_edges_x = np.histogram(x_coords, bins=256, range=(0, 255))
    bin_centers_x = (bin_edges_x[:-1] + bin_edges_x[1:]) / 2
    params_x, params_covariance_x = curve_fit(gauss, bin_centers_x, hist_x, p0=p0_x)
    errors_x = np.sqrt(np.diag(params_covariance_x))

    p0_y = [1000, np.mean(y_coords), 10]
    hist_y, bin_edges_y = np.histogram(y_coords, bins=256, range=(0, 255))
    bin_centers_y = (bin_edges_y[:-1] + bin_edges_y[1:]) / 2
    params_y, params_covariance_y = curve_fit(gauss, bin_centers_y, hist_y, p0=p0_y)
    errors_y = np.sqrt(np.diag(params_covariance_y))

    heatmap, xedges, yedges = np.histogram2d(x_coords, y_coords, bins=(256, 256), range=[[0, 255], [0, 255]])
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.cla()
    cax = ax.imshow(heatmap.T, origin='lower', cmap='viridis', norm=LogNorm())
    fig.colorbar(cax, ax=ax, label='Hits')
    ax.set_xlabel("Active pixels per event")
    ax.set_ylabel("Number of events")
    ax.grid(True)
    plt.savefig(direc + '/' + str(filename) + '.pdf', bbox_inches='tight', pad_inches=0.08)
    plt.close()

def main():
    # Get the arguments
    parser = argparse.ArgumentParser(description='Plot the occupancy of a run')
    parser.add_argument('runpath', type=str, help='Path to the hdf5 file or a folder containing multiple hdf5 files')
    parser.add_argument('--ecc', type=float, help='Selects a cut as a lower bound on the eccentricity of the event. The default is 1', default=1)
    args = parser.parse_args()

    run = args.runpath
    # Ploting if one file is provided
    if run.endswith('h5'):
        # Create the folder for storing the plot
        direc = os.path.dirname(run) + "/Occupancy"
        if os.path.exists(direc) == False:
            os.makedirs(direc)
        datafile = os.path.basename(run)
        f = h5py.File(run, 'r+')
        filename = datafile.replace('.h5', '')
        timepix_version = f['reconstruction'].attrs['TimepixVersion'][0].decode('utf-8')
        reconstruction = f['reconstruction']
        for name in reconstruction:
            x = f.get('reconstruction/' + name + '/chip_0/x')[:]
            y = f.get('reconstruction/' + name + '/chip_0/y')[:]
            x = np.concatenate(x)
            y = np.concatenate(y)
            ecc = f.get('reconstruction/' + name + '/chip_0/eccentricity')[:]
            if args.ecc > 1:
                hits = hits[ecc > args.ecc]
                filename = filename + '_ecc' + str(args.ecc)
            occupancy(x, y, direc, filename)
    # Plotting if a folder with files is provided
    elif os.path.isdir(run):
        print("Plotting the spectra for all hdf5 files in the folder")

        # Create the folder for storing the plots
        direc = os.path.dirname(run) + "/Occupancy"
        if os.path.exists(direc) == False:
            os.makedirs(direc)

        # Iterate over all hdf5 files in the folder
        files = [file for file in os.listdir(run) if file.endswith('.h5')]
        for file in tqdm(files):
            datei_pfad = os.path.join(run, file)
            with h5py.File(datei_pfad, 'r') as f:
                filename = file.replace('.h5', '')
                timepix_version = f['reconstruction'].attrs['TimepixVersion'][0].decode('utf-8')
                reconstruction = f['reconstruction']
                for name in reconstruction:
                    x = f.get('reconstruction/' + name + '/chip_0/x')[:]
                    y = f.get('reconstruction/' + name + '/chip_0/y')[:]
                    x = np.concatenate(x)
                    y = np.concatenate(y)
                    ecc = f.get('reconstruction/' + name + '/chip_0/eccentricity')[:]
                    if args.ecc > 1:
                        hits = hits[ecc > args.ecc]
                        filename = filename + '_ecc' + str(args.ecc)
                    occupancy(x, y, direc, filename)
    else:
        print("Please choose a correct data file or folder")    


if __name__ == "__main__":
    main()