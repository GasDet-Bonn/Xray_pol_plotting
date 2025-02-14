import numpy as np
from tqdm import tqdm
import os
import h5py
import argparse

import matplotlib.pyplot as plt

def calc_radius(coords, charges):
    x_pos, y_pos = coords
    center = np.average(coords, axis=1, weights=charges)
    X = np.vstack((x_pos - center[0], y_pos - center[1]))
    M = np.dot(X * charges, X.T)

    eigenvalues, eigenvectors = np.linalg.eig(M)
    principal_axis = eigenvectors[:, np.argmax(eigenvalues)]
    projection_xy_fit = np.dot(X[:2].T, principal_axis[:2])
    m2 = np.sum(charges * projection_xy_fit**2) / np.sum(charges)

    rad = np.sqrt((np.power(center[0]-127.273, 2) + np.power(center[1]-127.273, 2))/m2)
    return rad

# Create a histogram save it to the folder `direc` with the name `filename`.
def radius(rad, direc, filename):
    plt.clf()
    maximum = rad.max()
    plt.hist(rad, bins=100, range=(0, 10))
    plt.xlabel("Radius [pixel]")
    plt.ylabel("Number of events")
    plt.grid(True)
    plt.savefig(direc + '/' + str(filename) + '.pdf')
    #print(filename, np.mean(rad), np.std(rad))

def main():
    # Get the arguments
    parser = argparse.ArgumentParser(description='Plots a histogram of the distance of hits to the center of charge of an event')
    parser.add_argument('runpath', type=str, help='Path to the hdf5 file or a folder containing multiple hdf5 files')
    args = parser.parse_args()

    run = args.runpath
    # Ploting if one file is provided
    if run.endswith('h5'):
        # Create the folder for storing the plot
        direc = os.path.dirname(run) + "/Radius"
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
                x = f.get('reconstruction/' + name + '/chip_0/x')[:]
                y = f.get('reconstruction/' + name + '/chip_0/y')[:]
            except:
                print("No angles found - please perform the angular reconstruction")
                continue
            try:
                charge = f.get('reconstruction/' + name + '/chip_0/charge')[:]
            except:
                charge = f.get('reconstruction/' + name + '/chip_0/ToT')[:]
            rad = []
            for event in range(len(x)):
                rad.append(calc_radius(np.array([x[event], y[event]]), charge[event]))
            radius(np.array(rad), direc, filename)
    # Plotting if a folder with files is provided
    elif os.path.isdir(run):
        print("Plotting the spectra for all hdf5 files in the folder")

        # Create the folder for storing the plots
        direc = os.path.dirname(run) + "/Radius"
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
                        x = f.get('reconstruction/' + name + '/chip_0/x')[:]
                        y = f.get('reconstruction/' + name + '/chip_0/y')[:]
                    except:
                        print("No angles found - please perform the angular reconstruction for file", filename)
                        continue
                    try:
                        charge = f.get('reconstruction/' + name + '/chip_0/charge')[:]
                    except:
                        charge = f.get('reconstruction/' + name + '/chip_0/ToT')[:]
                    rad = []
                    for event in range(len(x)):
                        rad.append(calc_radius(np.array([x[event], y[event]]), charge[event]))
                    radius(np.array(rad), direc, filename)
    else:
        print("Please choose a correct data file or folder")    


if __name__ == "__main__":
    main()