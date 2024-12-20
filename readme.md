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
tqdm
h5py
```

## Scripts
The following scripts are currently part of the collection:

### Event plotting
Plots the active pixels per event together with the reconstructed angles. It can be
used as follows:
```
python3 plot_event.py [-h] [--pixel_size] [--png] [--events EVENTS] [--event_offset EVENT_OFFSET] [--plot3d] runpath
```
The `runpath` is a path to an Hdf5 file. Additionally there are the following optional parameters:
- `pixel_size` Use a variable pixel size for plotting based on the charge per pixel.
- `png` Save the plots as png in addition to the default pdf.
- `events`Select how many events are plotted starting with the first event. Default is 100 events. 0 for all events.
- `event_offset` Selects the first event that is plotted (Counting starts at 0). Default is 0.
- `plot3d` Add two additional plots in xz and yz projection. Only works with Timepix3 data.

### Pixel spectrum
Plots a histogram of the number of active pixels per event. It can be used as follows:
```
python3 plot_pixelspectrum.py [-h] runpath
```
The `runpath` can either be the full path to a hdf5 file but also to a folder containing
hdf5 files. In the latter case a plot for each hdf5 file is generated.

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

### Angles
Plots a histogram of reconstructed photo electron angles. It can be used as follows:
```
python3 plot_angles.py [-h] [--png] [--fit] [--log] runpath
```
The `runpath` can either be the full path to a hdf5 file but also to a folder containing
hdf5 files. The script only works if the files contain angular data. This data can be
reconstructed for example with [Xray_pol_reco](https://github.com/GasDet-Bonn/Xray_pol_reco).
There are the following optional parameters:
- `png` Save the plots as png in addition to the default pdf.
- `fit` A cos^2 fit is performed and added to the plot.
- `log` Write a log file with the modulation and the polarisation angle. Also activates the `fit`.
