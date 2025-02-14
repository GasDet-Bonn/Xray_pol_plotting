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

def lin(x, a, b):
    return a*x + b

# Create a histogram and save it to the folder `direc` with the name `filename`.
def chargeperpixel(charge, direc, filename, electrons):
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
    event_gain = np.array(list(map(np.nanmean, charge)))
    event_gain_error = np.array(list(map(np.nanstd, charge)))
    range_min, range_max = 0, 100

    ax.hist(event_gain, bins=101, range=(range_min, range_max), density=True)
    if electrons:
        ax.set_xlabel("Charge per pixel [electrons]")
    else:
        ax.set_xlabel("Charge per pixel [electrons]")
    ax.set_ylabel("Normalised number of events")
    ax.set_xlim(0, 100)
    #ax.set_yscale('log')
    plt.grid(True)
    plt.savefig(direc + '/' + str(filename) + '.pdf', bbox_inches='tight', pad_inches=0.08)
    plt.close()

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.cla()

    filter = event_gain_error > 0
    event_gain = event_gain[filter]
    event_gain_error = event_gain_error[filter]
    params, params_covariance = curve_fit(lin, range(len(event_gain)), event_gain, sigma=event_gain_error, p0=[0, 2000])
    params_errors = np.sqrt(np.diag(params_covariance))
    x_values = range(len(event_gain[::500]))
    lin_fit = lin(x_values, *params)

    # Calculate fit
    fit_values = lin(range(len(event_gain)), *params)

    # Calculate residuals
    residuals = event_gain - fit_values
    reduced_chi_squared = np.sum((residuals / event_gain_error) ** 2) / (len(event_gain) - len(params))

    fit_info = (f"$a = {params[0]:.3f} \\pm {params_errors[0]:.3f}$\n"
            f"$b = {params[1]:.3f} \\pm {params_errors[1]:.3f}$\n"
            f"$\\chi^2_{{\\mathrm{{red}}}} = {reduced_chi_squared:.3f}$")
    
    plt.gca().text(0.68, 0.95, fit_info, transform=plt.gca().transAxes, 
               fontsize=10, verticalalignment='top', horizontalalignment='left',
               bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.5))

    ax.errorbar(x=range(len(event_gain[::500])), y=event_gain[::500], yerr=event_gain_error[::500], fmt='o', zorder=1)
    ax.plot(x_values, lin_fit, 'r-', lw=2, zorder=2)
    if electrons:
        ax.set_ylabel("Mean charge per pixel [electrons]")
    else:
        ax.set_ylabel("Mean charge per pixel [electrons]")
    ax.set_xlabel("Event [100 events]")
    plt.grid(True)
    plt.savefig(direc + '/' + str(filename) + '_timespan.pdf', bbox_inches='tight', pad_inches=0.08)

def main():
    # Get the arguments
    parser = argparse.ArgumentParser(description='Plot the average charge per pixel per event. Only considers pixels of the start of the track for polarised data.')
    parser.add_argument('runpath', type=str, help='Path to the hdf5 file or a folder containing multiple hdf5 files')
    args = parser.parse_args()

    run = args.runpath
    electrons = False
    # Ploting if one file is provided
    if run.endswith('h5'):
        # Create the folder for storing the plot
        direc = os.path.dirname(run) + "/ChargeperpixelEventPol"
        if os.path.exists(direc) == False:
            os.makedirs(direc)
        datafile = os.path.basename(run)
        f = h5py.File(run, 'r+')
        filename = datafile.replace('.h5', '')
        timepix_version = f['reconstruction'].attrs['TimepixVersion'][0].decode('utf-8')
        reconstruction = f['reconstruction']
        for name in reconstruction:
            start_indices = f.get('reconstruction/' + name + '/chip_0/start_indices')[:]
            try:
                charge = f.get('reconstruction/' + name + '/chip_0/charge')[:]
                electrons = True
            except:
                charge = f.get('reconstruction/' + name + '/chip_0/ToT')[:]
            filtered_charge = [arr[indices.astype(int)] for arr, indices in zip(charge, start_indices)]
            filtered_charge = [arr for arr in filtered_charge if len(arr) > 0]
            chargeperpixel(filtered_charge, direc, filename, electrons)
    # Plotting if a folder with files is provided
    elif os.path.isdir(run):
        print("Plotting the spectra for all hdf5 files in the folder")

        # Create the folder for storing the plots
        direc = os.path.dirname(run) + "/ChargeperpixelEventPol"
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
                    start_indices = f.get('reconstruction/' + name + '/chip_0/start_indices')[:]
                    try:
                        charge = f.get('reconstruction/' + name + '/chip_0/charge')[:]
                        electrons = True
                    except:
                        charge = f.get('reconstruction/' + name + '/chip_0/ToT')[:]
                    filtered_charge = [arr[indices.astype(int)] for arr, indices in zip(charge, start_indices)]
                    filtered_charge = [arr for arr in filtered_charge if len(arr) > 1]
                    chargeperpixel(filtered_charge, direc, filename, electrons)
    else:
        print("Please choose a correct data file or folder")    


if __name__ == "__main__":
    main()