# Andor_kinetics
Created on Fri Mar 19 15:02:31 2021

@author: dd
version: 1.5

Extracts intensity decay from Andor .asc file series and saves results to 
kinetics.txt.

First, export Andor .sif data series as .asc files with parameters as follows: 
Tab>Signal, Separator>Comma, 
Append Aquisition Info>At top, Write each image as separate file

Each series must be saved in a seperate directory.

Integrated intensity is calculated from whole matrix and from region of interest.
ROI coordinates are in pixel numbers. ROI width x height must be set manually.
The center of ROI is selected in the first image with a single mouse click. 
Position of maximum intensity subarray of size width x height is displayed as 
a white rectangle in the first image.

Due to large number of files in the series, only user selected images will be 
displayed.
Use Sliders in the first figure to select file from the series.
Added two user slected images in kinetics window for visualizations.

Time delay between two images is read from info header as 
dt = cycle_time.

File numbering starts from '1'. 
Time delay counting starts from 'dt'.
 
'Background' intensity is set as (mean - 3*std.) from the last image and 
is subtracted from the whole series.

Intensity scale of images in the series is set to (mean +/- 2*std.dev.).

If ROI border exceeds edges of the matrix, new ROI center coordinates are 
retracted back into the matrix.

Issues:

Mouse input to work properly, change Spyder settings in
Preferences>IPhyton console>Graphics>Backend to Automatic.

ROI coordinates precision may vary by one point due to number 
rounding ("width/2").
"""
