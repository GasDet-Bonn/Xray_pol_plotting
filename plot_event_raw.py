import numpy as np
from tqdm import tqdm
import os
import h5py
import argparse

import matplotlib.pyplot as plt
import matplotlib.patches as patches

def rms(arr):
    return np.sqrt(np.mean(np.square(arr)))

def main():
    # Get the arguments
    parser = argparse.ArgumentParser(description='Plot raw event')
    parser.add_argument('--pixel_size', action='store_true', help='Use a variable pixel size for plotting based on the charge per pixel.')
    parser.add_argument('--png', action='store_true', help='Save the plots as png in addition to the default pdf.')
    parser.add_argument('--events', type=int, help='Select how many events are plotted starting with the first event. Default is 100 events. 0 for all events', default=100)
    parser.add_argument('--event_offset', type=int, help='Selects the first event that is plotted (Counting starts at 0). Default is 0.', default=0)
    parser.add_argument('runpath', type=str, help='Path to the hdf5 file')

    args = parser.parse_args()

    # Process the arguments
    run = args.runpath
    if run.endswith('h5'):
        datafile = os.path.basename(run)
    else:
        print("Please choose a correct data file")

    # Open the corresponding datafile
    f = h5py.File(run, 'r+')
    filename = datafile.replace('.h5', '')

    timepix_version = f['reconstruction'].attrs['TimepixVersion'][0].decode('utf-8')

    reconstruction = f['reconstruction']
    for name in reconstruction:
        # Load the x and y coordinates of the pixel with the option to rotate
        posx_raw = f.get('reconstruction/' + name + '/chip_0/x')[:]
        posy_raw = f.get('reconstruction/' + name + '/chip_0/y')[:]
        eventnumber = f.get('reconstruction/' + name + '/chip_0/eventNumber')[:]

        # If there is calibrated charge data use it, otherwise use the ToT
        try:
            charge = f.get('reconstruction/' + name + '/chip_0/charge')[:]
            chargetot = f.get('reconstruction/' + name + '/chip_0/ToT')[:]
        except:
            charge = f.get('reconstruction/' + name + '/chip_0/ToT')[:]
            chargetot = charge

        # Create a new folder based on the filename of the run
        direc = filename + "/Events_raw"
        if os.path.exists(direc) == False:
            os.makedirs(direc)

        # Set the range of events based on arguments
        if args.events == 0:
            eventrange = range(args.event_offset, len(eventnumber))
        else:
            if args.event_offset + args.events > len(eventnumber):
                raise ValueError(f"Error: event_offset+events exceeds the amount of available events {len(eventnumber)}.")
            eventrange = range(args.event_offset, args.event_offset + args.events)

        for event in tqdm(eventrange):
            try:
                fig, ax = plt.subplots()
                ax.set_xlabel('x [pixel]')
                ax.set_ylabel('y [pixel]')
                scatter = ax.scatter(x=posx_raw[event], y=posy_raw[event], c=charge[event], cmap='viridis', s=2, marker='s')
                colorbar = fig.colorbar(scatter, ax=ax, label='Charge [electrons]')
                ax.set_xlim(75, 175)
                ax.set_ylim(75, 175)

                # Set the aspect ratio for a quadratic plot
                ax.set_aspect('equal', adjustable='box')
                ax.grid(True)

                # Save the event as a pdf
                plt.savefig(direc + '/event-' + str(event) + '.pdf', bbox_inches='tight', pad_inches=0.08)
                # If selected save additionally as png
                if args.png:
                    plt.savefig(direc + '/event-' + str(event) + '.png', bbox_inches='tight', pad_inches=0.08)

            except Exception as e:
                print(f"Event {event} could not be plotted: {e}")


if __name__ == "__main__":
    main()