# -*- coding: utf-8 -*-
"""
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
    * mouse input to work properly, change Spyder settings in
    Preferences>IPhyton console>Graphics>Backend to Automatic.
    * ROI coordinates precision may vary by one point due to number 
    rounding ("width/2").
"""

import glob
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0})
from matplotlib.widgets import Slider, Button
from copy import deepcopy
import time


#path to series data, set manually
data_folder = 'C:\\DARBAI\\SOLIS\\Data\\'

#set roi's width and height manually
width = 30
height = 30

cmap = 'viridis' #<-- colormap for images
stdev = 2 #<-- number of standard dev. in intensity scale, set manually

int_total = []
int_roi = []
ttime = []


def findMaxSubmatrix(matrix, w, h):
    '''Finds maximum subarray of size w*h within a matrix (Kadane's algorithm)
        Parameters:
            w (int): An integer, width of subarray
            h (int): An integer, height of subarray
            matrix: 2d array
        Returns:
            best_pos (int): Top left corner of max subarray
            best (float): Maximum sum of subarray, #temporily removed.
    '''
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

def cnvrtInt(a):
    ''' Takes float, returns rounded and converted to int'''
    return int(round(a, 0))

def readNclean(filename):
    '''Prepares matrix (loads file, subtracts bg, applies CRR)'''
    dataframe = pd.read_csv(filename, skiprows=37, header=None)
    #remove first and last columns which are irrelevant
    dataframe.drop([0, len(dataframe.columns)-1], axis=1, inplace=True)
    # global background
    dataframe = dataframe - background
    #CRR
    dataframe[dataframe > 10*dataframe.values.flatten().std()] = dataframe.values.flatten().mean()
    return dataframe

#background 
last_filename = glob.glob(data_folder+'*.asc')[-1]
last = pd.read_csv(last_filename,skiprows=37, header=None)
last.drop([0, len(last.columns)-1], axis=1, inplace=True)
background = last.values.flatten().mean() - 3*last.values.flatten().std()

# plot the first image of the series to select region of interest (roi)
first_filename = glob.glob(data_folder+'*.asc')[0]
first = readNclean(first_filename)
best = findMaxSubmatrix(first.to_numpy(), width, height)

fig0, axs0 = plt.subplots(figsize=(6,6))
im0 = axs0.imshow(first, interpolation=None, cmap=cmap)

#plots maximum subarray of size WxH in a matrix
roi0 = plt.Rectangle(((best[1], best[0])),
                    width, height, fill=False, ec='white', linestyle='--')
plt.gca().add_patch(roi0)
axs0.scatter(best[1]+cnvrtInt(width/2), best[0]+cnvrtInt(height/2), s=2, c='k')

axs0.grid(True, linestyle='--', linewidth=0.5)
axs0.set_xlabel('Pixel number')
axs0.set_ylabel('Pixel number')
axs0.set_title('Center of max intensity {}x{} window is in:[{},{}]'.format(width,
                                                                    height,
                                                                    best[1]+cnvrtInt(width/2),
                                                                    best[0]+cnvrtInt(height/2)))
cbar = fig0.colorbar(im0, fraction=0.046, pad=0.04)
cbar.set_label('Intensity (a.u.)')

#set center of roi using mouse click
# mouse_in = plt.ginput(n=1)[0]
# centr = [cnvrtInt(x) for x in mouse_in]
centr = [15, 64]

tstart = time.time()
vcoord = centr[0]
hcoord = centr[1]

# border check
if hcoord-round(height/2)<0 or hcoord+round(height/2)>first.shape[1] or vcoord-round(width/2)<0 or hcoord+round(width/2)>first.shape[0]:
    print('Border detected!')
if hcoord-round(height/2)<0:
    hcoord=hcoord+abs(hcoord-round(height/2))
    print('New hcoord {}'.format(hcoord))
if hcoord+round(height/2)>first.shape[1]:
    hcoord=hcoord-(hcoord+round(height/2)-first.shape[1])
    print('New hcoord {}'.format(hcoord))
if vcoord-round(width/2)<0:
    vcoord=vcoord+abs(vcoord-round(width/2))
    print('New vcoord {}'.format(vcoord))
if vcoord+round(width/2)>first.shape[0]:
    vcoord=vcoord-(vcoord+round(width/2)-first.shape[1])
    print('New vcoord {}'.format(vcoord))        

#find time delay (dt) between two matrices
info=pd.read_csv(first_filename, nrows=12, header=None)
cycle_time = info.iloc[8].str.split(':').tolist()[0][1]
dt = float(cycle_time)
t=dt
n=0

#definitions for axes in subplots
left_ax, width_ax = 0.1, 0.65
bottom_ax, height_ax = 0.1, 0.65
spacing = 0.005
rect_image = [left_ax, bottom_ax, width_ax, height_ax]
rect_x = [left_ax, bottom_ax + height_ax + spacing, width_ax, 0.2]
rect_y = [left_ax + width_ax + spacing, bottom_ax, 0.2, height_ax]

#Slider to select file for intensity profile
axamp = plt.axes([0.25, 0.04, 0.50, 0.02])
samp = Slider(axamp, 'Intensity profile:', 0, len(glob.glob(data_folder+'*.asc')), 
              valinit=0, valstep=1)

def update_profile(val):
    '''Displays intensity profile in file selected by Slider'''
    # num is the current value of the slider padded with 0
    num = str(samp.val).zfill(4)
    # splits numbers from filename and adds back num from the slider
    file = glob.glob(data_folder+'*.asc')[0]
    file = data_folder + file.split('\\')[-1][:-8] + num + '.asc'

    #load the matrix
    df = readNclean(file)
    
    #plot the matrix
    fig1 = plt.figure(figsize=(6,6))
    axs1 = fig1.add_axes(rect_image)
    axs1.set_xlim(0, len(df.index))
    axs1.set_ylim(len(df.columns), 0)
    
    im1 = axs1.imshow(df, interpolation=None, cmap=cmap,
                vmin=df.values.flatten().mean()-stdev*df.values.flatten().std(),
                vmax=df.values.flatten().mean()+stdev*df.values.flatten().std())
    #user selected roi
    roi = plt.Rectangle((vcoord-width/2, hcoord-height/2),
                        width, height, fill=False, ec='red')
    plt.gca().add_patch(roi)
    axs1.axhline(y=hcoord, color='r', linestyle='--')
    axs1.axvline(x=vcoord, color='r', linestyle='--')
    axs1.text(len(df.index)+3, -24,
              '#:{} delay:{}s \ncenter:[{} {}] \n{}'.format(samp.val, 
                                                            round(samp.val*dt,1), 
                                                            vcoord, hcoord, 
                                                            file.split('\\')[-1]),
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

samp.on_changed(update_profile)

#Slider to select image
axamp2 = plt.axes([0.25, 0.01, 0.50, 0.02])
samp2 = Slider(axamp2, 'Single image:', 0, len(glob.glob(data_folder+'*.asc')), 
              valinit=0, valstep=1)

#Displays clean image of selected file
def update_img(val):
    '''Displays image without intensity profiles in files selected by Slider'''
    # num is the current value of the slider padded with 0
    num = str(samp2.val).zfill(4)
    # splits numbers from filename and adds back num from the slider
    file = glob.glob(data_folder+'*.asc')[0]
    file = data_folder + file.split('\\')[-1][:-8] + num + '.asc'
    #print(file)

    #load image
    img = readNclean(file)  
    
    fig, axs = plt.subplots(figsize=(6,6))
    im = axs.imshow(img, interpolation=None, cmap=cmap, 
                    vmin=img.values.flatten().mean()-stdev*img.values.flatten().std(), 
                    vmax=img.values.flatten().mean()+stdev*img.values.flatten().std())
    
    #user selected roi
    roi = plt.Rectangle((vcoord-width/2, hcoord-height/2),
                        width, height, fill=False, ec='red')
    plt.gca().add_patch(roi)
    
    axs.axis('off')
    axs.set_xlabel('Pixel number')
    axs.set_ylabel('Pixel number')
    axs.text(len(img.index)-25, -10, 
             '#:{} delay:{}s \n{}'.format(samp2.val,
                                          round(samp2.val*dt,1),
                                          file.split('\\')[-1]),
              bbox=dict(facecolor='blue', alpha=0.3))
    
    cbar = fig.colorbar(im, fraction=0.046, pad=0.04)
    cbar.set_label('Intensity (a.u.)')

samp2.on_changed(update_img)
    
#loop for the kinetics calculation
for file in glob.glob(data_folder+'*.asc'):
    #load the matrix
    df2 = readNclean(file)   
    
    #select region of interest in the matrix
    df3=df2.iloc[vcoord-cnvrtInt(width/2):vcoord+cnvrtInt(width/2), 
                hcoord-cnvrtInt(height/2):hcoord+cnvrtInt(height/2)]
    
    #calculate sum of the total matrix area and roi area
    intensity = df2.values.flatten().sum()
    int_total.append(intensity)
    intensity_roi = df3.values.flatten().sum()
    int_roi.append(intensity_roi)
    
    #time delta between matrices
    ttime.append(t)
    t+=dt
    n+=1

#plot the kinetics
fig2, axs2 = plt.subplots(figsize=(6,5))
axs2.scatter(ttime, int_total, label='Total area') #Normalized --> int_total/max(int_total
axs2.scatter(ttime, int_roi, label='ROI area') #Normalized --> int_roi/max(int_roi

axs2.set_xlim(left=-10)
axs2.tick_params(axis='both', which='both', direction='in')
axs2.set_xlabel('Delay (s)')
axs2.set_ylabel('Integrated Intensity (a.u.)')
axs2.legend(loc='best', frameon=False)
axs2.set_yscale('log')

#draws visualization images in kinetics window
#slider for image 1
axamp3 = fig2.add_axes([0.45, 0.97, 0.4, 0.02])
samp3 = Slider(axamp3, 'Image 1:', 0, len(glob.glob(data_folder+'*.asc')),
              valinit=0, valstep=1)
#slider for image 2
axamp4 = plt.axes([0.45, 0.92, 0.4, 0.02])
samp4 = Slider(axamp4, 'Image 2:', 0, len(glob.glob(data_folder+'*.asc')),
              valinit=0, valstep=1)
#clear button
rax = fig2.add_axes([0.13, 0.92, 0.17, 0.07])
check = Button(rax, ('Clear'))

#axis for visualizations
vizax1 = fig2.add_axes([0.4, 0.57, 0.25, 0.25], visible=False) 
vizax2 = fig2.add_axes([0.65, 0.57, 0.25, 0.25], visible=False) 

def update_img1(pic):
    '''Displays Image 1 in files selected by Slider
    Intensity scale of images is set to default for better contrast in 
    early stages
    '''
    num1 = str(samp3.val).zfill(4)
    # splits numbers from filename and adds back num from the slider
    file = glob.glob(data_folder+'*.asc')[0]
    file = data_folder + file.split('\\')[-1][:-8] + num1 + '.asc'
    #load image
    img1 = readNclean(file)
    vizax1.set_visible(True)
    vizax1.imshow(img1, cmap=cmap)
    vizax1.set_title('{}s'.format(round(int(num1)*dt, 1)))
    vizax1.axis('off')
    axs2.legend(loc='center right', frameon=False)
    plt.sca(vizax1)
    roi = plt.Rectangle((vcoord-width/2, hcoord-height/2),
                    width, height, fill=False, ec='red')
    plt.gca().add_patch(roi)    
    plt.draw()

samp3.on_changed(update_img1)

def update_img2(pic):
    '''Displays Image 2 in files selected by Slider
    Intensity scale of images is set to (mean +/- stdev*standard dev.) for 
    better contrast in later stages.
    '''
    num2 = str(samp4.val).zfill(4)
    # splits numbers from filename and adds back num from the slider
    file = glob.glob(data_folder+'*.asc')[0]
    file = data_folder + file.split('\\')[-1][:-8] + num2 + '.asc'
    #load image
    img2 = readNclean(file)
    vizax2.set_visible(True)
    vizax2.imshow(img2, cmap=cmap,
                  vmin=img2.values.flatten().mean()-stdev*img2.values.flatten().std(), 
                  vmax=img2.values.flatten().mean()+stdev*img2.values.flatten().std())
    vizax2.set_title('{}s'.format(round(int(num2)*dt, 1)))
    vizax2.axis('off')
    axs2.legend(loc='center right', frameon=False)
    plt.sca(vizax2)
    roi = plt.Rectangle((vcoord-width/2, hcoord-height/2),
                width, height, fill=False, ec='red')
    plt.gca().add_patch(roi)     
    plt.draw()
    
samp4.on_changed(update_img2)

def clear(self):
    '''Clears Image 1 and Image 2'''
    axs2.legend(loc='best', frameon=False)
    try:
        vizax1.cla()
    except NameError:
        pass
    try:
        vizax2.cla()
    except NameError:
        pass
    try:
        vizax1.axis('off')
    except NameError:
        pass
    try:
        vizax2.axis('off')
    except NameError:
        pass
check.on_clicked(clear)

#histogram of the last image
# fig3, ax3 = plt.subplots()
# ax3.hist(df)

#save to file
meta = ('roi_center: {} {} \nwidth: {} \nheight: {}'.format(vcoord, hcoord, 
                                                        width, height))
result = pd.DataFrame(list(zip(ttime, int_total, int_roi)), 
                            columns =['delay', 'total_int', 'roi_int'])
result.to_csv(data_folder+'kinetics.txt', index=False)
with open(data_folder+'kinetics.txt', 'a') as f:
    f.write(meta)
    
print('Elapsed time: {} s'.format(round(time.time()-tstart,2)))