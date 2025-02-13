import numpy as np
from tqdm import tqdm
import os
import h5py
import argparse

import matplotlib.pyplot as plt
import matplotlib.patches as patches

def main():
    # Get the arguments
    parser = argparse.ArgumentParser(description='Plot events with reconstructed angles')
    parser.add_argument('--pixel_size', action='store_true', help='Use a variable pixel size for plotting based on the charge per pixel.')
    parser.add_argument('--png', action='store_true', help='Save the plots as png in addition to the default pdf.')
    parser.add_argument('--events', type=int, help='Select how many events are plotted starting with the first event. Default is 100 events. 0 for all events', default=100)
    parser.add_argument('--event_offset', type=int, help='Selects the first event that is plotted (Counting starts at 0). Default is 0.', default=0)
    parser.add_argument('--plot3d', action='store_true', help='Add two additional plots in xz and yz projection. Only works with Timepix3 data.')
    parser.add_argument('--radius', type=float, nargs=2, help='Plot with radius')
    parser.add_argument('--weight', type=float, help='Plot with weighting')
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

    if timepix_version == 'Timepix1' and args.plot3d:
        print('The `plot3d` option works only with Timepix3 data - ignoring it.')

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
            theta1 = f.get('reconstruction/' + name + '/chip_0/theta_firststage')[:]
            theta2 = f.get('reconstruction/' + name + '/chip_0/theta_secondstage')[:]

        if args.radius:
            absorption_point_x = f.get('reconstruction/' + name + '/chip_0/absorption_point_x')[:]
            absorption_point_y = f.get('reconstruction/' + name + '/chip_0/absorption_point_y')[:]

        # If there is calibrated charge data use it, otherwise use the ToT
        try:
            charge = f.get('reconstruction/' + name + '/chip_0/charge')[:]
        except:
            charge = f.get('reconstruction/' + name + '/chip_0/ToT')[:]

        if args.plot3d and timepix_version == 'Timepix3':
            fig, (ax, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
        else:
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
                charges_event = charge[event]
                center_x_all = np.average(posx_raw[event], weights=charges_event)
                center_y_all = np.average(posy_raw[event], weights=charges_event)
                center_x_start = np.average(posx_raw[event][start_indices[event].astype(int)], weights=charges_event[start_indices[event].astype(int)])
                center_y_start = np.average(posy_raw[event][start_indices[event].astype(int)], weights=charges_event[start_indices[event].astype(int)])

                line_length = 20
                ax.cla()
                ax.set_xlabel('x [pixel]')
                ax.set_ylabel('y [pixel]')
                
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
                if not args.radius:
                    if args.pixel_size:
                        max_charge = np.max(charges_event)
                        ax.scatter(posx_raw[event][start_indices[event].astype(int)], posy_raw[event][start_indices[event].astype(int)], color='orange', marker='s', s=4*charge[event][start_indices[event].astype(int)]/max_charge)
                        ax.scatter(posx_raw[event][end_indices[event].astype(int)], posy_raw[event][end_indices[event].astype(int)], color='blue', marker='s', s=4*charge[event][end_indices[event].astype(int)]/max_charge)
                    else:
                        ax.scatter(posx_raw[event][start_indices[event].astype(int)], posy_raw[event][start_indices[event].astype(int)], color='orange', marker='s', s=2)
                        ax.scatter(posx_raw[event][end_indices[event].astype(int)], posy_raw[event][end_indices[event].astype(int)], color='blue', marker='s', s=2)
                else:
                    if args.pixel_size:
                        max_charge = np.max(charges_event)
                        ax.scatter(posx_raw[event][start_indices[event].astype(int)], posy_raw[event][start_indices[event].astype(int)], color='orange', marker='s', s=4*charge[event][start_indices[event].astype(int)]/max_charge)
                        ax.scatter(posx_raw[event][end_indices[event].astype(int)], posy_raw[event][end_indices[event].astype(int)], color='orange', marker='s', s=4*charge[event][end_indices[event].astype(int)]/max_charge)
                    else:
                        ax.scatter(posx_raw[event][start_indices[event].astype(int)], posy_raw[event][start_indices[event].astype(int)], color='orange', marker='s', s=2)
                        ax.scatter(posx_raw[event][end_indices[event].astype(int)], posy_raw[event][end_indices[event].astype(int)], color='orange', marker='s', s=2)

                # Plot the center of charge of all pixels
                ax.scatter(center_x_all, center_y_all, color='red', label="Center of charge", marker='o', s=40)

                if not args.radius:
                    # Plot the center of charge for the start pixels
                    ax.scatter(center_x_start, center_y_start, color='green', label="Center of charge", marker='o', s=40)
                else:
                    ax.scatter(absorption_point_x[event], absorption_point_y[event], color='green', label="Absorption point", marker='o', s=40)

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
                
                if args.radius:
                    try:
                        posz_raw = ((toa[event] * 25 + 1) - ftoa[event] * 1.5625)*3.59/55.
                    except:
                        posz_raw = np.zeros(len(posx_raw[event]))
                    coords = np.array([posx_raw[event], posy_raw[event], posz_raw])
                    center = np.average(coords, axis=1, weights=charges_event)
                    X = np.vstack((posx_raw[event] - center[0], posy_raw[event] - center[1], posz_raw - center[2]))
                    M = np.dot(X * charges_event, X.T)
                    eigenvalues, eigenvectors = np.linalg.eig(M)
                    principal_axis = eigenvectors[:, np.argmax(eigenvalues)]
                    projection_new_plane = np.dot(X.T, principal_axis)
                    m2_max = np.sum(charges_event*np.power(projection_new_plane, 2)) / np.sum(charges_event)
                    radius_min, radius_max = args.radius
                    radius_low = radius_min * np.sqrt(m2_max)
                    radius_high = radius_max * np.sqrt(m2_max)
                    circle_low = plt.Circle((center_x_all, center_y_all), radius_low, color='r', fill=False, linestyle='--')
                    circle_high = plt.Circle((center_x_all, center_y_all), radius_high, color='r', fill=False, linestyle='--')
                    ax.add_patch(circle_low)
                    ax.add_patch(circle_high)

                if not args.radius:
                # Create an arrow for the reconstructed angle with only the start pixels
                    x1 = center_x_start - line_length*np.cos(phi2[event])
                    x2 = center_x_start + line_length*np.cos(phi2[event])
                    y1 = center_y_start - line_length*np.sin(phi2[event])
                    y2 = center_y_start + line_length*np.sin(phi2[event])
                else:
                    x1 = absorption_point_x[event] - line_length*np.cos(phi2[event])
                    x2 = absorption_point_x[event] + line_length*np.cos(phi2[event])
                    y1 = absorption_point_y[event] - line_length*np.sin(phi2[event])
                    y2 = absorption_point_y[event] + line_length*np.sin(phi2[event])
                arrow = patches.FancyArrowPatch((x1, y1), (x2, y2), mutation_scale=5, arrowstyle='->,head_width=1,head_length=2', color='green')
                ax.add_patch(arrow)

                # Set the aspect ratio for a quadratic plot
                ax.set_aspect('equal', adjustable='box')
                ax.grid(True)

                if args.plot3d and timepix_version == 'Timepix3':
                    ax2.cla()
                    ax3.cla()
                    posz_raw = ((toa[event] * 25 + 1) - ftoa[event] * 1.5625)#*velocity/55.

                    center_z_all = np.average(posz_raw, weights=charges_event)
                    center_z_start = np.average(posz_raw[start_indices[event].astype(int)], weights=charges_event[start_indices[event].astype(int)])
                    # Plot the center of charge of all pixels
                    ax2.scatter(center_x_all, center_z_all, color='red', label="Center of charge", marker='o', s=40)

                    # Plot the center of charge for the start pixels
                    ax2.scatter(center_x_start, center_z_start, color='green', label="Center of charge", marker='o', s=40)
                    ax2.scatter(posx_raw[event][start_indices[event].astype(int)], posz_raw[start_indices[event].astype(int)], color='black', marker='s', s=2)
                    ax2.scatter(posx_raw[event][end_indices[event].astype(int)], posz_raw[end_indices[event].astype(int)], color='blue', marker='s', s=2)

                    ax2.set_xlabel('x [pixel]')
                    ax2.set_ylabel('z [ns]')
                    ax2.grid(True)

                    # Crate a dashed orhogonal line to the arrow based on all pixels
                    x1 = center_x_all + line_length*np.sin(theta1[event])
                    x2 = center_x_all - line_length*np.sin(theta1[event])
                    y1 = center_z_all - line_length*np.cos(theta1[event])
                    y2 = center_z_all + line_length*np.cos(theta1[event])
                    ax2.plot([x1, x2], [y1, y2], linestyle='--', color='red')

                    ax3.scatter(center_y_all, center_z_all, color='red', label="Center of charge", marker='o', s=40)

                    # Plot the center of charge for the start pixels
                    ax3.scatter(center_y_start, center_z_start, color='green', label="Center of charge", marker='o', s=40)
                    ax3.scatter(posy_raw[event][start_indices[event].astype(int)], posz_raw[start_indices[event].astype(int)], color='black', marker='s', s=2)
                    ax3.scatter(posy_raw[event][end_indices[event].astype(int)], posz_raw[end_indices[event].astype(int)], color='blue', marker='s', s=2)

                    ax3.set_xlabel('y [pixel]')
                    ax3.set_ylabel('z [ns]')
                    ax3.grid(True)

                    # Set ranges for the plots. For x and y the same range as for the xy plot is used while for z a new is defined
                    if np.max(posx_raw[event]) - np.min(posx_raw[event]) > np.max(posy_raw[event]) - np.min(posy_raw[event]):
                        ax2.set_xlim(np.min(posx_raw[event]) - 20, np.max(posx_raw[event]) + 20)
                        x_range = np.max(posx_raw[event]) - np.min(posx_raw[event])
                        ax3.set_xlim(center_y_all - x_range/2-20, center_y_all + x_range/2+20)
                    else:
                        ax3.set_xlim(np.min(posy_raw[event]) - 20, np.max(posy_raw[event]) + 20)
                        y_range = np.max(posy_raw[event]) - np.min(posy_raw[event])
                        ax2.set_xlim(center_x_all - y_range/2-20, center_x_all + y_range/2+20)
                    zrange = np.max(posz_raw) - np.min(posz_raw)
                    ax2.set_ylim(np.min(posz_raw)-zrange/2, np.max(posz_raw)+zrange/2)
                    ax3.set_ylim(np.min(posz_raw)-zrange/2, np.max(posz_raw)+zrange/2)

                    plt.tight_layout()

                # Save the event as a pdf
                plt.savefig(direc + '/event-' + str(event) + '.pdf')
                # If selected save additionally as png
                if args.png:
                    plt.savefig(direc + '/event-' + str(event) + '.png')

            except Exception as e:
                print(f"Event {event} could not be plotted: {e}")


if __name__ == "__main__":
    main()