import numpy as np
from tqdm import tqdm
import os
import h5py
import argparse

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def gauss(x, mu, sigma, A):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

def rms(arr):
    return np.sqrt(np.mean(np.square(arr)))

# Create a histogram and save it to the folder `direc` with the name `filename`.
def totrms(charge, direc, filename, fit=False):
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

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    ax.cla()
    chargerms = np.array(list(map(rms, charge)))
    if fit:
        counts, bins, _ = plt.hist(chargerms, bins=100, range=(0, 300))
        bin_centers = (bins[:-1] + bins[1:]) / 2
        popt, pcov = curve_fit(gauss, bin_centers, counts, p0=[150, 10, 1], bounds=([-np.inf, 0, -np.inf], [np.inf, np.inf, np.inf]), maxfev=100000)
        mu, sigma, A = popt
        perr = np.sqrt(np.diag(pcov))

        fit_info = (f"$\\mu = ({mu:.3f} \\pm {perr[0]:.3f})\ cc$\n"
                f"$\\sigma = ({sigma:.3f} \\pm {perr[1]:.3f})\ cc$\n"
                f"$a = {A:.0f} \\pm {perr[2]:.0f}$")
        
        plt.gca().text(0.68, 0.95, fit_info, transform=plt.gca().transAxes, 
                fontsize=10, verticalalignment='top', horizontalalignment='left',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.5))
    
        x = np.linspace(0, 300, 1000)
        plt.plot(x, gauss(x, *popt), 'k', linewidth=2)
        print(filename, mu-2*sigma, mu+2*sigma)
    else:
        ax.hist(chargerms, bins=100, range=(0, 300))
    #plt.hist(charge, bins=int(maximum*1.1)+1, range=(0, maximum*1.1))
    ax.set_xlabel("ToT RMS [clock cycles]")
    #ax.set_xlim(0, 250000)
    ax.set_ylabel("Number of events")
    plt.grid(True)
    plt.savefig(direc + '/' + str(filename) + '.pdf', bbox_inches='tight', pad_inches=0.08)

def main():
    # Get the arguments
    parser = argparse.ArgumentParser(description='Plot the ToT rms')
    parser.add_argument('runpath', type=str, help='Path to the hdf5 file or a folder containing multiple hdf5 files')
    args = parser.parse_args()

    run = args.runpath
    electrons = False
    # Ploting if one file is provided
    if run.endswith('h5'):
        # Create the folder for storing the plot
        direc = os.path.dirname(run) + "/TOTrms"
        if os.path.exists(direc) == False:
            os.makedirs(direc)
        datafile = os.path.basename(run)
        f = h5py.File(run, 'r+')
        filename = datafile.replace('.h5', '')
        timepix_version = f['reconstruction'].attrs['TimepixVersion'][0].decode('utf-8')
        reconstruction = f['reconstruction']
        for name in reconstruction:
            charge = f.get('reconstruction/' + name + '/chip_0/ToT')[:]
            totrms(charge, direc, filename, True)
    # Plotting if a folder with files is provided
    elif os.path.isdir(run):
        print("Plotting the spectra for all hdf5 files in the folder")

        # Create the folder for storing the plots
        direc = os.path.dirname(run) + "/TOTrms"
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
                    charge = f.get('reconstruction/' + name + '/chip_0/ToT')[:]
                    totrms(charge, direc, filename, True)
    else:
        print("Please choose a correct data file or folder")    


if __name__ == "__main__":
    main()