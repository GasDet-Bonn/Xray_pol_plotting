# Collection of plotting scripts for GridPix-based X-ray polarimetry

This repository is a collection of various plotting scripts for X-ray polarimetry with
GridPix detectors based on Timepix and Timepix3. By default the script require Hdf5 files
that are already processed by [TimepixAnalysis](https://github.com/Vindaar/TimepixAnalysis)
and the polarimetry [reconstruction](https://github.com/GasDet-Bonn/Xray_pol_reco). If
requirements per script deviate from this it is stated individually.

## Requirements
All scripts require the following python packages:
```
numpy
scipy
tqdm
h5py
```

## Scripts
The following scripts are currently part of the collection:

### Raw Event plotting
Plots the active pixels per event without polarisation angles. It can be
used as follows:
```
python3 plot_event_raw.py [-h] [--pixel_size] [--png] [--events EVENTS] [--event_offset EVENT_OFFSET] runpath
```
The `runpath` is a path to an Hdf5 file. Additionally there are the following optional parameters:
- `pixel_size` Use a variable pixel size for plotting based on the charge per pixel.
- `png` Save the plots as png in addition to the default pdf.
- `events`Select how many events are plotted starting with the first event. Default is 100 events. 0 for all events.
- `event_offset` Selects the first event that is plotted (Counting starts at 0). Default is 0.


### Event plotting
Plots the active pixels per event together with the reconstructed angles. It can be
used as follows:
```
python3 plot_event.py [-h] [--pixel_size] [--png] [--events EVENTS] [--event_offset EVENT_OFFSET] [--plot3d] [--radius INNER OUTER] [--weight WEIGHT] runpath
```
The `runpath` is a path to an Hdf5 file. Additionally there are the following optional parameters:
- `pixel_size` Use a variable pixel size for plotting based on the charge per pixel.
- `png` Save the plots as png in addition to the default pdf.
- `events`Select how many events are plotted starting with the first event. Default is 100 events. 0 for all events.
- `event_offset` Selects the first event that is plotted (Counting starts at 0). Default is 0.
- `plot3d` Add two additional plots in xz and yz projection. Only works with Timepix3 data.
- `radius` Add circles for the inner and outer radii of the weight based reconstructions
- `weight` Provide the weight for the weight based reconstruction


### Pixel spectrum
Plots a histogram of the number of active pixels per event. It can be used as follows:
```
python3 plot_pixelspectrum.py [-h] [--ecc] runpath
```
The `runpath` can either be the full path to a hdf5 file but also to a folder containing
hdf5 files. In the latter case a plot for each hdf5 file is generated.
There is the following optional parameter:
- `ecc` Selects a cut as a lower bound on the eccentricity of the event. The default is 1.


### Charge spectrum
Plots a histogram of the total charge per event. It can be used as follows:
```
python3 plot_chargespectrum.py [-h] runpath
```
The `runpath` can either be the full path to a hdf5 file but also to a folder containing
hdf5 files. In the latter case a plot for each hdf5 file is generated. If a hdf5 file
contains calibrated charge data in electrons this is used, otherwise the ToT is used.


### Charge per pixel
Plots a histogram of the charge per pixel per event. It can be used as follows:
```
python3 plot_chargeperpixel.py [-h] runpath
```
The `runpath` can either be the full path to a hdf5 file but also to a folder containing
hdf5 files. In the latter case a plot for each hdf5 file is generated. If a hdf5 file
contains calibrated charge data in electrons this is used, otherwise the ToT is used.


### Average charge per pixel
Plots the average charge per pixel per event and fits the result with a linear function.
It can be used as follows:
```
python3 plot_chargeperpixel_event.py [-h] runpath
```
The `runpath` can either be the full path to a hdf5 file but also to a folder containing
hdf5 files. In the latter case a plot for each hdf5 file is generated. If a hdf5 file
contains calibrated charge data in electrons this is used, otherwise the ToT is used.


### Average charge per pixel - track start
Plots the average charge per pixel per event and fits the result with a linear function.
Only the hits that are considered as part of the start of a photo electron track are
considered. Thus, a dataset with a polarisation reconstruction is required.
It can be used as follows:
```
python3 plot_chargeperpixel_event.py [-h] runpath
```
The `runpath` can either be the full path to a hdf5 file but also to a folder containing
hdf5 files. In the latter case a plot for each hdf5 file is generated. If a hdf5 file
contains calibrated charge data in electrons this is used, otherwise the ToT is used.

### Angles
Plots a histogram of reconstructed photo electron angles. It can be used as follows:
```
python3 plot_angles.py [-h] [--png] [--fit] [--log] [--ecc] runpath
```
The `runpath` can either be the full path to a hdf5 file but also to a folder containing
hdf5 files. The script only works if the files contain angular data. This data can be
reconstructed for example with [Xray_pol_reco](https://github.com/GasDet-Bonn/Xray_pol_reco).
There are the following optional parameters:
- `png` Save the plots as png in addition to the default pdf.
- `fit` A cos^2 fit is performed and added to the plot.
- `log` Write a log file with the modulation and the polarisation angle. Also activates the `fit`.
- `ecc` Selects a cut as a lower bound on the eccentricity of the event. The default is 1.


### Occupancy
Plots a 2D histogram of hits per pixel of a full run. It can be used as follows:
```
python3 plot_occupancy.py [-h] [--ecc] runpath
```
The `runpath` can either be the full path to a hdf5 file but also to a folder containing
hdf5 files. In the latter case a plot for each hdf5 file is generated. 
There is the following optional parameter:
- `ecc` Selects a cut as a lower bound on the eccentricity of the event. The default is 1.


### Pixel - Total charge
Plots a correlation of the number of hits in an event and the total collected charge of an event.
It can be used as follows:
```
python3 plot_pixel_totalcharge.py [-h] runpath
```
The `runpath` can either be the full path to a hdf5 file but also to a folder containing
hdf5 files. In the latter case a plot for each hdf5 file is generated. 


### Charge per pixel - Total charge
Plots a correlation of the the total charge in an event and the charges of the individual hits in
the event. It can be used as follows:
```
python3 plot_chargeperpixel_totalcharge.py [-h] runpath
```
The `runpath` can either be the full path to a hdf5 file but also to a folder containing
hdf5 files. In the latter case a plot for each hdf5 file is generated. 


### ToT RMS
Plots the RMS of the ToT in events. It can be used as follows:
```
python3 plot_totrms.py [-h] runpath
```
The `runpath` can either be the full path to a hdf5 file but also to a folder containing
hdf5 files. In the latter case a plot for each hdf5 file is generated.


### Radius
Plots the distance of hits to the center of charge of their respective event. It can be used as follows:
```
python3 plot_radius.py [-h] runpath
```
The `runpath` can either be the full path to a hdf5 file but also to a folder containing
hdf5 files. In the latter case a plot for each hdf5 file is generated.
