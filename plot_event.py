import numpy as np
from tqdm import tqdm
import os
import h5py
import argparse

import matplotlib.pyplot as plt
import matplotlib.patches as patches

def main():
    # Get the arguments
    parser = argparse.ArgumentParser(description='Plot reconstructed polarisation angles per event')
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
        phi1 = f.get('reconstruction/' + name + '/chip_0/angle_firststage')[:]
        phi2 = f.get('reconstruction/' + name + '/chip_0/angle_secondstage')[:]
        eventnumber = f.get('reconstruction/' + name + '/chip_0/eventNumber')[:]
        start_indices = f.get('reconstruction/' + name + '/chip_0/start_indices')[:]
        end_indices = f.get('reconstruction/' + name + '/chip_0/end_indices')[:]

        # Fot Timepix3 use additionally the timestamp per pixel
        if timepix_version == 'Timepix3':
            toa = f.get('reconstruction/' + name + '/chip_0/ToA')[:]
            ftoa = f.get('reconstruction/' + name + '/chip_0/fToA')[:]

        # If there is calibrated charge data use it, otherwise use the ToT
        try:
            charge = f.get('reconstruction/' + name + '/chip_0/charge')[:]
        except:
            charge = f.get('reconstruction/' + name + '/chip_0/ToT')[:]

        fig, ax = plt.subplots()

        # Create a new folder based on the filename of the run
        direc = filename + "/Events"
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
                # Calculate charge centers
                center_x_all = np.average(posx_raw[event], weights=charge[event])
                center_y_all = np.average(posy_raw[event], weights=charge[event])
                center_x_start = np.average(posx_raw[event][start_indices[event].astype(int)], weights=charge[event][start_indices[event].astype(int)])
                center_y_start = np.average(posy_raw[event][start_indices[event].astype(int)], weights=charge[event][start_indices[event].astype(int)])

                line_length = 20
                ax.cla()
                plt.xlabel('x [pixel]')
                plt.ylabel('y [pixel]')
                
                # Calculate the limits of the axis, results in a quadratic range
                if np.max(posx_raw[event]) - np.min(posx_raw[event]) > np.max(posy_raw[event]) - np.min(posy_raw[event]):
                    ax.set_xlim(np.min(posx_raw[event]) - 20, np.max(posx_raw[event]) + 20)
                    x_range = np.max(posx_raw[event]) - np.min(posx_raw[event])
                    ax.set_ylim(center_y_all - x_range/2-20, center_y_all + x_range/2+20)
                else:
                    ax.set_ylim(np.min(posy_raw[event]) - 20, np.max(posy_raw[event]) + 20)
                    y_range = np.max(posy_raw[event]) - np.min(posy_raw[event])
                    ax.set_xlim(center_x_all - y_range/2-20, center_x_all + y_range/2+20)

                # Plot the start and end pixels. Either with a static pixel size or a size based on the normalised charge per pixel
                if args.pixel_size:
                    max_charge = np.max(charge[event])
                    ax.scatter(posx_raw[event][start_indices[event].astype(int)], posy_raw[event][start_indices[event].astype(int)], color='black', marker='s', s=4*charge[event][start_indices[event].astype(int)]/max_charge)
                    ax.scatter(posx_raw[event][end_indices[event].astype(int)], posy_raw[event][end_indices[event].astype(int)], color='blue', marker='s', s=4*charge[event][end_indices[event].astype(int)]/max_charge)
                else:
                    ax.scatter(posx_raw[event][start_indices[event].astype(int)], posy_raw[event][start_indices[event].astype(int)], color='black', marker='s', s=2)
                    ax.scatter(posx_raw[event][end_indices[event].astype(int)], posy_raw[event][end_indices[event].astype(int)], color='blue', marker='s', s=2)

                # Plot the center of charge of all pixels
                ax.scatter(center_x_all, center_y_all, color='red', label="Center of charge", marker='o', s=40)

                # Plot the center of charge for the start pixels
                ax.scatter(center_x_start, center_y_start, color='green', label="Center of charge", marker='o', s=40)

                # Crate an arrow for the reconstructed angle with all pixels
                x1 = center_x_all - line_length*np.cos(phi1[event])
                x2 = center_x_all + line_length*np.cos(phi1[event])
                y1 = center_y_all - line_length*np.sin(phi1[event])
                y2 = center_y_all + line_length*np.sin(phi1[event])
                arrow = patches.FancyArrowPatch((x1, y1), (x2, y2), mutation_scale=5, arrowstyle='->,head_width=1,head_length=2', color='red')
                ax.add_patch(arrow)

                # Crate a dashed orhogonal line to the arrow based on all pixels
                x1 = center_x_all + line_length*np.sin(phi1[event])
                x2 = center_x_all - line_length*np.sin(phi1[event])
                y1 = center_y_all - line_length*np.cos(phi1[event])
                y2 = center_y_all + line_length*np.cos(phi1[event])
                ax.plot([x1, x2], [y1, y2], linestyle='--', color='red')

                # Create an arrow for the reconstructed angle with only the start pixels
                x1 = center_x_start - line_length*np.cos(phi2[event])
                x2 = center_x_start + line_length*np.cos(phi2[event])
                y1 = center_y_start - line_length*np.sin(phi2[event])
                y2 = center_y_start + line_length*np.sin(phi2[event])
                arrow = patches.FancyArrowPatch((x1, y1), (x2, y2), mutation_scale=5, arrowstyle='->,head_width=1,head_length=2', color='green')
                ax.add_patch(arrow)

                # Set the aspect ratio for a quadratic plot
                ax.set_aspect('equal', adjustable='box')
                ax.grid(True)

                # Save the event as a pdf
                plt.savefig(direc + '/event-' + str(event) + '.pdf')
                # If selected save additionally as png
                if args.png:
                    plt.savefig(direc + '/event-' + str(event) + '.png')

            except Exception as e:
                print(f"Event {event} could not be plotted: {e}")


if __name__ == "__main__":
    main()