import numpy as np
from tqdm import tqdm
import os
import h5py
import argparse

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def cos_squared(x, A, B, C):
    return A * np.cos(x - C)**2 + B

# Create a histogram based on a array of angles and save it
# to the folder `direc` with the name `filename`.
def angles(angle, direc, filename, png, fit, log):
    plt.clf()

    if fit or log:
        angle = angle[~np.isnan(angle)]
        counts, bin_edges = np.histogram(angle, bins=100)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        errors = np.sqrt(counts)

        params, covariance = curve_fit(cos_squared, bin_centers, counts, bounds=([0, 0, -np.pi], [np.inf, np.inf, np.pi]), maxfev=1000000)
        A, B, C = params

        expected_counts = cos_squared(bin_centers, *params)
        chi_square = np.sum(((counts - expected_counts) ** 2) / errors ** 2)
        dof = len(counts) - len(params)
        reduced_chi_square = chi_square / dof

        errors = np.sqrt(np.diag(covariance))
        A_error, B_error, C_error = errors
        modulationfactor = A / (A + 2*B)
        modulationfactor_error = modulationfactor * np.sqrt((A_error / A)**2 + (2 * B_error / (A + 2 * B))**2)
        modulationfactor = modulationfactor * 100
        modulationfactor_error = modulationfactor_error * 100

    maximum = angle.max()
    plt.hist(angle, bins=100, range=(-np.pi, np.pi))
    if fit or log:
        plt.plot(bin_centers, cos_squared(bin_centers, A, B, C), label='Fit: A*cos^2(x) + B', color='red')
        ax = plt.gca()
        default_fontsize = ax.xaxis.get_label().get_size()
        text = f'$\mu = ({modulationfactor:.2f} \pm {modulationfactor_error:.2f}) \%$\n$\phi = ({C:.2f} \pm {C_error:.2f})$'
        if log:
            with open(direc + '/modulations.txt', 'a') as file:
                file.write(filename + '\t' + str(modulationfactor) + '\t' + str(modulationfactor_error) + '\t' + str(C) + '\t' + str(C_error) + '\t' + str(efficiency) + '\t' + str(reduced_chi_square) + '\t' + str(len(angle)) + '\n')

    plt.xlabel("Reconstructed angle [rad]")
    plt.ylabel("Number of events")
    plt.xlim(-np.pi, np.pi)
    plt.grid(True)
    plt.savefig(direc + '/' + str(filename) + '.pdf')
    if png:
        plt.savefig(direc + '/' + str(filename) + '.png')

def main():
    # Get the arguments
    parser = argparse.ArgumentParser(description='Plot the pixelspectra (active pixels per event)')
    parser.add_argument('runpath', type=str, help='Path to the hdf5 file or a folder containing multiple hdf5 files')
    parser.add_argument('--png', action='store_true', help='Save the plots as png in addition to the default pdf.')
    parser.add_argument('--fit', action='store_true', help='Perform a cos^2 fit to the data.')
    parser.add_argument('--log', action='store_true', help='Write modulation factor and angle to a log file. Also activates the fit.')
    args = parser.parse_args()

    run = args.runpath
    # Ploting if one file is provided
    if run.endswith('h5'):
        # Create the folder for storing the plot
        direc = os.path.dirname(run) + "/Angle"
        if os.path.exists(direc) == False:
            os.makedirs(direc)
        datafile = os.path.basename(run)
        f = h5py.File(run, 'r+')
        filename = datafile.replace('.h5', '')
        timepix_version = f['reconstruction'].attrs['TimepixVersion'][0].decode('utf-8')
        reconstruction = f['reconstruction']
        for name in reconstruction:
            try:
                angle = f.get('reconstruction/' + name + '/chip_0/angle_secondstage')[:]
            except:
                print("No angles found - please perform the angular reconstruction")
            angles(angle, direc, filename, args.png, args.fit, args.log)
    # Plotting if a folder with files is provided
    elif os.path.isdir(run):
        print("Plotting the spectra for all hdf5 files in the folder")

        # Create the folder for storing the plots
        direc = os.path.dirname(run) + "/Angle"
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
                        angle = f.get('reconstruction/' + name + '/chip_0/angle_secondstage')[:]
                    except:
                        print("No angles found - please perform the angular reconstruction for file", filename)
                    angles(angle, direc, filename, args.png, args.fit, args.log)
    else:
        print("Please choose a correct data file or folder")    


if __name__ == "__main__":
    main()