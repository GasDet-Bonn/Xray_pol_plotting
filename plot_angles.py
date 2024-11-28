import numpy as np
from tqdm import tqdm
import os
import h5py
import argparse

import matplotlib.pyplot as plt

# Create a histogram based on a array of angles and save it
# to the folder `direc` with the name `filename`.
def angles(angle, direc, filename):
    plt.clf()
    maximum = angle.max()
    plt.hist(angle, bins=100, range=(-np.pi, np.pi))
    plt.xlabel("Reconstructed angle [rad]")
    plt.ylabel("Number of events")
    plt.xlim(-np.pi, np.pi)
    plt.grid(True)
    plt.savefig(direc + '/' + str(filename) + '.pdf')

def main():
    # Get the arguments
    parser = argparse.ArgumentParser(description='Plot the pixelspectra (active pixels per event)')
    parser.add_argument('runpath', type=str, help='Path to the hdf5 file or a folder containing multiple hdf5 files')
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
            angles(angle, direc, filename)
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
                    angles(angle, direc, filename)
    else:
        print("Please choose a correct data file or folder")    


if __name__ == "__main__":
    main()