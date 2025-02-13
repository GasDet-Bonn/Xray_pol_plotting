import numpy as np
from tqdm import tqdm
import os
import h5py
import argparse
import scipy
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt

def polya(x, K, G, T):
    Tp1 = T + 1
    return (K / G) * (np.power(Tp1, Tp1))/(scipy.special.gamma(Tp1)) * np.power((x/G), T) * np.exp(-(Tp1*x)/(G))

# Create a histogram based on a array of charges per pixel and save it
# to the folder `direc` with the name `filename`.
def chargeperpixel(charge, direc, filename, electrons):
    # Definition of the plot size
    fig_width, fig_height = 7, 6  # in inches

    # Relative font size
    font_size = fig_height * 2  # Beispiel: 2 mal die Breite der Figur

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
    maximum = charge.max()
    unique_data = np.sort(np.unique(charge))
    bin_edges = np.concatenate((unique_data - 0.5, [unique_data[-1] + 0.5]))

    hist, bins = np.histogram(charge, bins=bin_edges, density=True)
    bin_centers = (bins[:-1] + bins[1:]) / 2

    fit_range_start = 600
    fit_range = (bin_centers >= fit_range_start)

    params, params_covariance = curve_fit(polya, bin_centers[fit_range], hist[fit_range], p0=[10000, 2000, 1])
    params_errors = np.sqrt(np.diag(params_covariance))
    x_values = np.linspace(fit_range_start, max(bin_edges), 1000)
    polya_fit = polya(x_values, *params)

    residuals = hist[fit_range] - polya(bin_centers[fit_range], *params)
    reduced_chi_squared = np.sum((residuals ** 2) / polya(bin_centers[fit_range], *params)) / (len(hist[fit_range]) - len(params))

    fit_info = (f"$K = {params[0]:.3f} \\pm {params_errors[0]:.3f}$\n"
            f"$G = {params[1]:.3f} \\pm {params_errors[1]:.3f}$\n"
            f"$\\theta = {params[2]:.3f} \\pm {params_errors[2]:.3f}$\n"
            f"$\\chi^2_{{\\mathrm{{red}}}} = {reduced_chi_squared:.3f}$")
    
    plt.gca().text(0.68, 0.95, fit_info, transform=plt.gca().transAxes, 
               fontsize=10, verticalalignment='top', horizontalalignment='left',
               bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.5))

    ax.hist(charge, bins=bin_edges, density=True)
    plt.plot(x_values, polya_fit, 'r-', lw=2)
    if electrons:
        ax.set_xlabel("Charge per pixel [electrons]")
    else:
        ax.set_xlabel("Charge per pixel [clock cycles]")
    ax.set_ylabel("Normalised number of events")
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
        direc = os.path.dirname(run) + "/Chargeperpixel"
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
                electrons = True
            except:
                charge = f.get('reconstruction/' + name + '/chip_0/ToT')[:]
            charge = np.concatenate(charge)
            chargeperpixel(charge, direc, filename, electrons)
    # Plotting if a folder with files is provided
    elif os.path.isdir(run):
        print("Plotting the spectra for all hdf5 files in the folder")

        # Create the folder for storing the plots
        direc = os.path.dirname(run) + "/Chargeperpixel"
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
                        electrons = True
                    except:
                        charge = f.get('reconstruction/' + name + '/chip_0/ToT')[:]
                    charge = np.concatenate(charge)
                    chargeperpixel(charge, direc, filename, electrons)
    else:
        print("Please choose a correct data file or folder")    


if __name__ == "__main__":
    main()