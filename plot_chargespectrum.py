import numpy as np
from tqdm import tqdm
import os
import h5py
import argparse

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def gauss(x, mu, sigma, A):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

# Create a histogram based on a array of total charge per event and save it
# to the folder `direc` with the name `filename`.
def chargespectrum(charge, charge_tot, direc, filename, electrons, fit=False):
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
    if fit:
        counts, bins, _ = plt.hist(charge, bins=100, range=(10000, 500000), alpha=1.0)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        popt, pcov = curve_fit(gauss, bin_centers, counts, p0=[100000, 5000, 50], bounds=([-np.inf, 0, -np.inf], [np.inf, np.inf, np.inf]), maxfev=100000)
        mu, sigma, A = popt
        perr = np.sqrt(np.diag(pcov))
        x = np.linspace(1000, 500000, 5000)
        ax.plot(x, gauss(x, *popt), 'k', linewidth=2)
        fwhm = 2 * np.sqrt(2 * np.log(2)) * sigma
        fwhm_error = 2 * np.sqrt(2 * np.log(2)) * perr[1]
        res = 100*sigma/mu
        res_error = 100 * np.sqrt(np.power(perr[1]/mu, 2) + np.power(perr[0]*sigma/(mu*mu), 2))
        res_fwhm = 100*fwhm/mu
        res_fwhm_error = 100 * np.sqrt(np.power(fwhm_error/mu, 2) + np.power(perr[0]*fwhm/(mu*mu), 2))
        print(filename, mu, perr[0], sigma, perr[1], fwhm, fwhm_error, res, res_error, res_fwhm, res_fwhm_error)

        fit_info = (f"$\\mu = ({mu:.0f} \\pm {perr[0]:.0f})\ e^-$\n"
                f"$\\sigma = ({sigma:.0f} \\pm {perr[1]:.0f})\ e^-$\n"
                f"$a = {A:.0f} \\pm {perr[2]:.0f}$")
        
        plt.gca().text(0.68, 0.95, fit_info, transform=plt.gca().transAxes, 
                fontsize=10, verticalalignment='top', horizontalalignment='left',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.5))
    else:
        ax.hist(charge, bins=100, range=(0, 1000000))
    if electrons:
        ax.set_xlabel("Total charge per event [electrons]")
    else:
        ax.set_xlabel("Total charge per event [clock cycles]")
    ax.set_ylabel("Number of events")
    plt.grid(True)
    plt.savefig(direc + '/' + str(filename) + '.pdf', bbox_inches='tight', pad_inches=0.08)

def main():
    # Get the arguments
    parser = argparse.ArgumentParser(description='Plot the pixelspectra (active pixels per event)')
    parser.add_argument('runpath', type=str, help='Path to the hdf5 file or a folder containing multiple hdf5 files')
    args = parser.parse_args()

    run = args.runpath
    electrons = False
    # Ploting if one file is provided
    if run.endswith('h5'):
        # Create the folder for storing the plot
        direc = os.path.dirname(run) + "/Chargespectrum"
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
            charge = np.array(list(map(np.sum, charge)))
            chargespectrum(charge, charge_tot, direc, filename, electrons, True)
    # Plotting if a folder with files is provided
    elif os.path.isdir(run):
        print("Plotting the spectra for all hdf5 files in the folder")

        # Create the folder for storing the plots
        direc = os.path.dirname(run) + "/Chargespectrum"
        if os.path.exists(direc) == False:
            os.makedirs(direc)

        # Iterate over all hdf5 files in the folder
        files = [file for file in os.listdir(run) if file.endswith('.h5')]
        for file in tqdm(files):
            datei_pfad = os.path.join(run, file)
            try:
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
                        charge = np.array(list(map(np.sum, charge)))
                    chargespectrum(charge, charge_tot, direc, filename, electrons, True)
            except:
                print(file, 'error')
    else:
        print("Please choose a correct data file or folder")    


if __name__ == "__main__":
    main()