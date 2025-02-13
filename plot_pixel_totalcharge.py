import numpy as np
from tqdm import tqdm
import os
import h5py
import argparse
import scipy
from scipy.optimize import curve_fit
import matplotlib.colors as mcolors

import matplotlib.pyplot as plt

def polya(x, K, G, T):
    Tp1 = T + 1
    return (K / G) * (np.power(Tp1, Tp1))/(scipy.special.gamma(Tp1)) * np.power((x/G), T) * np.exp(-(Tp1*x)/(G))

def rms(arr):
    return np.sqrt(np.mean(np.square(arr)))

# Create a histogram and save it to the folder `direc` with the name `filename`.
def pixels_totalcharge(charge, charge_tot, direc, filename, electrons):
    # Definition of the plot size
    fig_width, fig_height = 7, 6  # in inches

    chargetotal = np.array(list(map(np.sum, charge)))
    hits = np.array(list(map(len, charge)))

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

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    ax.cla()

    plt.hist2d(hits, chargetotal, bins=[100, 100], range=[[0,300],[0,1000000]], cmap='viridis', norm=mcolors.LogNorm())
    plt.colorbar(label='Number of hits')
    ax.set_xlabel("Hits per event [electrons]")
    ax.set_ylabel("Total charge per event [electrons]")
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 1000000)
    plt.grid(True)
    plt.savefig(direc + '/' + str(filename) + '.pdf', bbox_inches='tight', pad_inches=0.08)

def main():
    # Get the arguments
    parser = argparse.ArgumentParser(description='Plot a correlation of hits and total charge per event')
    parser.add_argument('runpath', type=str, help='Path to the hdf5 file or a folder containing multiple hdf5 files')
    args = parser.parse_args()

    run = args.runpath
    electrons = False
    # Ploting if one file is provided
    if run.endswith('h5'):
        # Create the folder for storing the plot
        direc = os.path.dirname(run) + "/PixelTotalcharge"
        if os.path.exists(direc) == False:
            os.makedirs(direc)
        datafile = os.path.basename(run)
        f = h5py.File(run, 'r+')
        filename = datafile.replace('.h5', '')
        timepix_version = f['reconstruction'].attrs['TimepixVersion'][0].decode('utf-8')
        reconstruction = f['reconstruction']
        for name in reconstruction:
            try:
                charge = f.get('reconstruction/' + name + '/chip_0/charge')[:]
                charge_tot = f.get('reconstruction/' + name + '/chip_0/ToT')[:]
                electrons = True
            except:
                charge = f.get('reconstruction/' + name + '/chip_0/ToT')[:]
                charge_tot = charge
            pixels_totalcharge(charge, charge_tot, direc, filename, electrons)
    # Plotting if a folder with files is provided
    elif os.path.isdir(run):
        print("Plotting the spectra for all hdf5 files in the folder")

        # Create the folder for storing the plots
        direc = os.path.dirname(run) + "/PixelTotalcharge"
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
                    try:
                        charge = f.get('reconstruction/' + name + '/chip_0/charge')[:]
                        charge_tot = f.get('reconstruction/' + name + '/chip_0/ToT')[:]
                        electrons = True
                    except:
                        charge = f.get('reconstruction/' + name + '/chip_0/ToT')[:]
                        charge_tot = charge
                    pixels_totalcharge(charge, charge_tot, direc, filename, electrons)
    else:
        print("Please choose a correct data file or folder")    


if __name__ == "__main__":
    main()