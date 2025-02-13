import numpy as np
from tqdm import tqdm
import os
import h5py
import argparse

import matplotlib.pyplot as plt

# Create a histogram based on a array of hits per event and save it
# to the folder `direc` with the name `filename`.
def pixelspectrum(hits, direc, filename):
    # Definition der Abbildungsgröße
    fig_width, fig_height = 7, 6  # in inches

    # Anpassung der Schriftgröße relativ zur Abbildungsgröße
    font_size = fig_height * 2  # Beispiel: 2 mal die Breite der Figur

    # Einstellung von rcParams für konsistente Schriftgröße
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
    #maximum = hits.max()
    plt.clf()
    maximum = hits.max()
    ax.hist(hits, bins=int(maximum*1.1)+1, range=(0, maximum*1.1))
    ax.set_xlabel("Active pixels per event")
    ax.set_ylabel("Number of events")
    ax.grid(True)
    plt.savefig(direc + '/' + str(filename) + '.pdf', bbox_inches='tight', pad_inches=0.08)

def main():
    # Get the arguments
    parser = argparse.ArgumentParser(description='Plot the pixelspectra (active pixels per event)')
    parser.add_argument('runpath', type=str, help='Path to the hdf5 file or a folder containing multiple hdf5 files')
    parser.add_argument('--ecc', type=float, help='Selects a cut as a lower bound on the eccentricity of the event. The default is 1', default=1)
    args = parser.parse_args()

    run = args.runpath
    # Ploting if one file is provided
    if run.endswith('h5'):
        # Create the folder for storing the plot
        direc = os.path.dirname(run) + "/Pixelspectrum"
        if os.path.exists(direc) == False:
            os.makedirs(direc)
        datafile = os.path.basename(run)
        f = h5py.File(run, 'r+')
        filename = datafile.replace('.h5', '')
        timepix_version = f['reconstruction'].attrs['TimepixVersion'][0].decode('utf-8')
        reconstruction = f['reconstruction']
        for name in reconstruction:
            hits = f.get('reconstruction/' + name + '/chip_0/hits')[:]
            ecc = f.get('reconstruction/' + name + '/chip_0/eccentricity')[:]
            if args.ecc > 1:
                hits = hits[ecc > args.ecc]
                filename = filename + '_ecc' + str(args.ecc)
            pixelspectrum(hits, direc, filename)
    # Plotting if a folder with files is provided
    elif os.path.isdir(run):
        print("Plotting the spectra for all hdf5 files in the folder")

        # Create the folder for storing the plots
        direc = os.path.dirname(run) + "/Pixelspectrum"
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
                    hits = f.get('runs/' + name + '/chip_0/Hits')[:]
                    ecc = f.get('reconstruction/' + name + '/chip_0/eccentricity')[:]
                    if args.ecc > 1:
                        hits = hits[ecc > args.ecc]
                        filename = filename + '_ecc' + str(args.ecc)
                    pixelspectrum(hits, direc, filename)
    else:
        print("Please choose a correct data file or folder")    


if __name__ == "__main__":
    main()