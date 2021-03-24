# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 15:02:31 2021

@author: dd
version: 1.1

Extracts intensity decay from Andor .asc files and saves results to kinetics.txt.

First export Andor .sif data series as .asc files with parameters as follows: 
Tab>Signal, Separator>Comma, 
Append Aquisition Info>At top, Write each image as separate file

Each series must be saved in a seperate directory.

Integrated intensity is calculated from whole matrix and from region of interest.
ROI coordinates are in pixel numbers. ROI width x height must be set manually.
The center of ROI is selected in the first image with a single mouse click. 
Position of maximum intensity subarray of size width x height is displayed as 
a white rectangle in the first image.

Time delay between two images is calculated automatically as 
di = exposure_time * number_of_accumulations.
 
'Background' intensity is set as minimum pixels' value from the last image and 
is subtracted from the whole series.

Intensity scale of images in the series is set to (mean +/- 2*std.dev.).

Issues:
    * mouse input to work properly, change Spyder settings in
    Preferences>IPhyton console>Graphics>Backend to Automatic.
    * ROI coordinates presicion may vary by one point due to number rounding ("width/2").
"""

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0})
from copy import deepcopy

#path to series data, set manually
data_folder = 'C:\\DARBAI\\SOLIS\\Data\\'

#set roi's width and height manually
width = 15
height = 30

cmap = 'viridis' #<-- colormap for images
stdev = 2 #<-- number of standard dev. in intensity scale, set manually

int_total = []
int_roi = []
delta = []

#finds maximum subarray of size width x height within a matrix (Kadane's algorithm)
def findMaxSubmatrix(matrix, w, h):
    nrows = len(matrix)
    ncols = len(matrix[0])           

    cumulative_sum = deepcopy(matrix)

    for r in range(nrows):
        for c in range(ncols):
            if r == 0 and c == 0:
                cumulative_sum[r][c] = matrix[r][c]
            elif r == 0:
                cumulative_sum[r][c] = cumulative_sum[r][c-1] + matrix[r][c]
            elif c == 0:
                cumulative_sum[r][c] = cumulative_sum[r-1][c] + matrix[r][c]
            else:
                cumulative_sum[r][c] = cumulative_sum[r-1][c] + cumulative_sum[r][c-1] - cumulative_sum[r-1][c-1] + matrix[r][c]

    best = 0
    best_pos = None

    for r1 in range(nrows):
        for c1 in range(ncols):
            r2 = r1 + h - 1
            c2 = c1 + w - 1
            if r2 >= nrows or c2 >= ncols:
                continue
            if r1 == 0 and c1 == 0:
                sub_sum = cumulative_sum[r2][c2]
            elif r1 == 0:
                sub_sum = cumulative_sum[r2][c2] - cumulative_sum[r2][c1-1]
            elif c1 == 0:
                sub_sum = cumulative_sum[r2][c2] - cumulative_sum[r1-1][c2]
            else:
                sub_sum = cumulative_sum[r2][c2] - cumulative_sum[r1-1][c2] - cumulative_sum[r2][c1-1] + cumulative_sum[r1-1][c1-1]
            if best < sub_sum:
                best_pos = r1,c1
                best = sub_sum

    # print ("maximum sum is:", best)
    # print ("top left corner on:", best_pos)
    return best_pos

#rounds float and converts to integer
def cnvrtInt(a):
    return int(round(a, 0))

# plot the first image of the series to select region of interest (roi)
first_filename = glob.glob(data_folder+'*.asc')[0]
first = pd.read_csv(first_filename, skiprows=37, header=None)
first.drop([0, len(first.columns)-1], axis=1, inplace=True)
best = findMaxSubmatrix(first.to_numpy(), width, height)

fig0, axs0 = plt.subplots()
im0 = axs0.imshow(first, interpolation=None, cmap=cmap)
#plots maximum subarray of size WxH in a matrix
roi0 = plt.Rectangle(((best[1], best[0])),
                    width, height, fill=False, ec='white', linestyle='--')
plt.gca().add_patch(roi0)
axs0.scatter(best[1]+cnvrtInt(width/2), best[0]+cnvrtInt(height/2), s=2, c='k')

axs0.grid(True, linestyle='--', linewidth=0.5)
axs0.set_xlabel('Pixel number')
axs0.set_ylabel('Pixel number')
axs0.set_title('Max intensity in '+str(width)+'x'+str(height)+' window is at:'+str(best[1]+cnvrtInt(width/2))+','+str(best[0]+cnvrtInt(height/2)))
cbar = fig0.colorbar(im0)
cbar.set_label('Intensity (a.u.)')

#set center of roi using mouse click
mouse_in = plt.ginput(n=1)[0]
centr = [cnvrtInt(x) for x in mouse_in]
# centr = (40, 65)

vcoord = centr[0]
hcoord = centr[1]

#calculate time delay (di) between two matrices as exposure_time*num_of_acc
info=pd.read_csv(first_filename, nrows=12, header=None)
exposure_time = info.iloc[7].str.split(':').tolist()[0][1]
num_of_acc = info.iloc[11].str.split(':').tolist()[0][1]
di = float(exposure_time)*float(num_of_acc)
i=0
n=0 

#finds min value in last image for background subtraction
last_filename = glob.glob(data_folder+'*.asc')[-1]
last = pd.read_csv(last_filename,skiprows=37, header=None)
last.drop([0, len(last.columns)-1], axis=1, inplace=True)
background = min(last.values.flatten())
# background = 4955

#definitions for axes in subplots
left_ax, width_ax = 0.1, 0.65
bottom_ax, height_ax = 0.1, 0.65
spacing = 0.005
rect_image = [left_ax, bottom_ax, width_ax, height_ax]
rect_x = [left_ax, bottom_ax + height_ax + spacing, width_ax, 0.2]
rect_y = [left_ax + width_ax + spacing, bottom_ax, 0.2, height_ax]

for file in os.listdir(data_folder):
    if file.endswith(".asc"):
        path = os.path.join(data_folder, file)
        
        #load the matrix
        df = pd.read_csv(path, skiprows=37, header=None)
        
        #remove first and last columns which are irrelevant
        df.drop([0, len(df.columns)-1], axis=1, inplace=True)
        
        #subtract background value
        df = df - background
        
        #plot all matrices
        fig1 = plt.figure(figsize=(6,6))
        axs1 = fig1.add_axes(rect_image)
        axs1.set_xlim(0, len(df.index))
        axs1.set_ylim(len(df.columns), 0)
        im1 = axs1.imshow(df, interpolation=None, cmap=cmap,
                    vmin=df.values.flatten().mean()-stdev*df.values.flatten().std(),
                    vmax=df.values.flatten().mean()+stdev*df.values.flatten().std())
        roi = plt.Rectangle((centr[0]-width/2, centr[1]-height/2),
                            width, height, fill=False, ec='red')
        plt.gca().add_patch(roi)
        axs1.axhline(y=hcoord, color='r', linestyle='--')
        axs1.axvline(x=vcoord, color='r', linestyle='--')
        axs1.text(len(df.index)+3, -24,
                  '#:'+str(n)+' delay:'+str(round(i,1))+'s'+
                  '\n'+'center:'+str(centr)+'\n'+file,
                  bbox=dict(facecolor='blue', alpha=0.3))
        
        ax_x = fig1.add_axes(rect_x, sharex=axs1)
        ax_x.tick_params(axis="x", labelbottom=False)
        #plot single line of horizontal crossection
        ax_x.plot(df.iloc[hcoord], color='#1f77b4', label='horiz')
        ax_x.legend(loc='lower left', fontsize='x-small')
        #plot integrated horizontal crossection
        ax_x2 = ax_x.twinx()
        ax_x2.plot(df.iloc[hcoord-int(height/2):hcoord+int(height/2)].sum(axis=0),
                  color='#ff7f0e', label=str(height)+'pts horiz')
        ax_x2.axis('off')
        ax_x2.legend(loc='lower right', fontsize='x-small')
        ax_x.axvline(x=vcoord-int(width/2), color='r')
        ax_x.axvline(x=vcoord+int(width/2), color='r')

        ax_y = fig1.add_axes(rect_y, sharey=axs1)
        ax_y.tick_params(axis="y", labelleft=False)
        #plot single line of vertical crossection
        ax_y.plot(df.iloc[:,vcoord], df.columns, color='#1f77b4', label='vert')
        ax_y.legend(loc='upper right', fontsize='x-small')
        #plot integrated vertical crossection
        ax_y2 = ax_y.twiny()
        ax_y2.plot(df.iloc[:,vcoord-int(width/2):vcoord+int(width/2)].sum(axis=1),
                  df.columns, color='#ff7f0e',label=str(width)+'pts vert')
        ax_y2.axis('off')
        ax_y2.legend(loc='lower right', fontsize='x-small')        
        ax_y.axhline(y=hcoord-int(height/2), color='r')
        ax_y.axhline(y=hcoord+int(height/2), color='r')
        
        #select region of interest in the matrix
        df2=df.iloc[centr[0]-cnvrtInt(width/2):centr[0]+cnvrtInt(width/2), 
                    centr[1]-cnvrtInt(height/2):centr[1]+cnvrtInt(height/2)]
        
        #calculate sum of the total matrix and roi
        intensity = df.values.flatten().sum()
        int_total.append(intensity)
        intensity_roi = df2.values.flatten().sum()
        int_roi.append(intensity_roi)
        
        #time delta between matrices
        delta.append(i)
        i+=di
        n+=1

#plot the kinetics
fig2, axs2 = plt.subplots()
axs2.scatter(delta, int_total/max(int_total), label='full image')
axs2.scatter(delta, int_roi/max(int_roi), label='roi image')
axs2.set_xlabel('Delay (s)')
axs2.set_ylabel('Normalized Intensity (a.u.)')
axs2.legend(loc='best')
axs2.set_yscale('log')
axs2.grid(True)

# #save to file
meta = ('center:'+str(centr)+'\n'+'width:'+str(width)+'\n'+'height:'+str(height))
result = pd.DataFrame(list(zip(delta, int_total, int_roi)), 
                            columns =['delay', 'total_int', 'roi_int'])
result.to_csv(data_folder+'kinetics.txt')
with open(data_folder+'kinetics.txt', 'a') as f:
    f.write(meta)
