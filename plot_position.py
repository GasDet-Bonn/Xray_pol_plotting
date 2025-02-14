import numpy as np
from tqdm import tqdm
import os
import h5py
import argparse

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def gauss(x, mu, sigma, A):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

def calc_center(coords, charges):
    x_pos, y_pos = coords
    coords = coords - 127.273
    center = np.average(coords, axis=1, weights=charges)
    distance = np.sqrt(np.power(center[0], 2) + np.power(center[1], 2))*55
    return distance

def rms(arr):
    return np.sqrt(np.mean(np.square(arr)))

# Create a histogram and save it to the folder `direc` with the name `filename`.
def position(pos, direc, filename, charge, charge_tot, fit=False):
    plt.clf()

    if fit:
        counts, bins, _ = plt.hist(pos, bins=100, range=(-30, 30), density=True, alpha=1)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        popt, pcov = curve_fit(gauss, bin_centers, counts, p0=[0, 10, 1], bounds=([-np.inf, 0, -np.inf], [np.inf, np.inf, np.inf]))
        mu, sigma, A = popt
        perr = np.sqrt(np.diag(pcov))
        x = np.linspace(-30, 30, 1000)
        plt.plot(x, gauss(x, *popt), 'k', linewidth=2)
    else:
        plt.hist(pos, bins=100, range=(-30, 30))
    plt.xlabel("Position [Pixel]")
    plt.ylabel("Number of events")
    plt.grid(True)
    plt.savefig(direc + '/' + str(filename) + '.pdf')

    if fit:
        return mu, sigma, A, perr
    else:
        return

def main():
    # Get the arguments
    parser = argparse.ArgumentParser(description='Plot the distance of the reconstructed absorption point to the beamspot')
    parser.add_argument('runpath', type=str, help='Path to the hdf5 file or a folder containing multiple hdf5 files')
    parser.add_argument('--beamspot', nargs=2, type=float, help='x and y coordinates of the beamspot')
    args = parser.parse_args()
    beamspot_x, beamspot_y = args.beamspot

    run = args.runpath
    # Ploting if one file is provided
    if run.endswith('h5'):
        # Create the folder for storing the plot
        direc = os.path.dirname(run) + "/Position"
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
                x = f.get('reconstruction/' + name + '/chip_0/absorption_point_x')[:]
                y = f.get('reconstruction/' + name + '/chip_0/absorption_point_y')[:]
            except:
                print("No angles found - please perform the angular reconstruction")
                continue

            distance = np.sqrt(np.power(x-beamspot_x, 2) + np.power(y-beamspot_y, 2))
            position(distance, direc, filename)
            mean = np.nanmean(distance)
            error = np.nanstd(distance)
            print(filename, mean, error)

    # Plotting if a folder with files is provided
    elif os.path.isdir(run):
        print("Plotting the spectra for all hdf5 files in the folder")

        # Create the folder for storing the plots
        direc = os.path.dirname(run) + "/Position"
        if os.path.exists(direc) == False:
            os.makedirs(direc)

        # Iterate over all hdf5 files in the folder
        files = [file for file in os.listdir(run) if file.endswith('.h5')]
        for file in files:
            datei_pfad = os.path.join(run, file)
            with h5py.File(datei_pfad, 'r') as f:
                filename = file.replace('.h5', '')
                timepix_version = f['reconstruction'].attrs['TimepixVersion'][0].decode('utf-8')
                reconstruction = f['reconstruction']
                for name in reconstruction:
                    try:
                        angle = f.get('reconstruction/' + name + '/chip_0/angle_secondstage')[:]
                        x = f.get('reconstruction/' + name + '/chip_0/absorption_point_x')[:]
                        y = f.get('reconstruction/' + name + '/chip_0/absorption_point_y')[:]
                    except:
                        print("No angles found - please perform the angular reconstruction for file", filename)
                        continue

                    distance = np.sqrt(np.power(x-beamspot_x, 2) + np.power(y-beamspot_y, 2))
                    position(distance, direc, filename)
                    mean = np.nanmean(distance)
                    error = np.nanstd(distance)

                    with open(direc + '/distance.txt', 'a') as file:
                        file.write(filename + '\t' + str(mean) + '\t' + str(error) + '\n')
    else:
        print("Please choose a correct data file or folder")    


if __name__ == "__main__":
    main()