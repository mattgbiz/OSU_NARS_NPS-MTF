#MTF Calculation written by Matt Bisbee 2/3/21
#This program takes the Edge Spread Function from Imagej as an input then calculates the LSF and MTF
#I save the input data as a CSV
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft,fftfreq
import cv2,os
import scipy as scipy
from scipy import optimize
from matplotlib import gridspec
from scipy import signal

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def CalculateMTF(ImagePath,EffectivePixelSize,ROI_Input,RoiNeeded=True):
    sensorLp = 1/(2*EffectivePixelSize) #lp/mm for sensor
    PMAG = 1    #assuming primary magnification size is 1
    lpPermm = sensorLp*PMAG         #cut-off frequency nyquist
    PixelPermm = 2*lpPermm       #number of pixels per mm
    #should PixelPermm just be 1/effective pixel size? so basically 10? what is it now?
    print('Current pixels per mm is {}'.format(PixelPermm))
    InputImg = cv2.imread(ImagePath,-1)
    #We are going to take an ROI from the image, need it to be brightened to see
    BitMax = 65536 #not actually the bit maximum but is close
    Imgmax = np.max(InputImg)
    print(Imgmax)
    StretchedImg = (InputImg*int(BitMax/Imgmax))
    print('select ROI for Edge'' left edge')
    if RoiNeeded == True:
        ImgROI = cv2.selectROI('Select Edge',StretchedImg)  #roi output is x,y coord of top left most point, then width is how far to right from that, height how far downc
        print(ImgROI)
    else:
        ImgROI = ROI_Input
    #if the user types c, this means cancel so we can just put in "expected" ROI 
    if not np.any((ImgROI)):
        print('contains only zeros, using defaults')
        ImgROI = np.array([219,188,56,162])

    #the order here is weird as we need rows then columns meaning y_start to y_start+height then x_start to x_start+width
    ImgCrop = InputImg[ImgROI[1]:(ImgROI[1]+ImgROI[3]),ImgROI[0]:(ImgROI[0]+ImgROI[2])]
    print(ImgCrop)
    #Now get the ESF from the cropped data
    Esf = np.mean(ImgCrop, axis=0)
    xvals = np.arange(1,len(Esf)+1,1)
    #normalize the ESF
    Ymax = np.max(Esf)
    normESF = Esf/Ymax

    #calculate the derivative of ESF to get LSF
    Lsf = np.diff(normESF)/np.diff(xvals)*-1
    #But we lose an element so need an X2 that is between each of these
    i = 0
    dxvals = np.zeros(len(Lsf))
    for i in range(len(Lsf)):
        dxvals[i] = (xvals[i+1]+xvals[i])/2

    #now get the mtf
    #number of sampling points
    N = len(Lsf)
    #Sample spacing
    Fs = PixelPermm  #sampling frequency
    Leng = len(xvals)   #how many pixels we have
    #need the next power of 2, aka figure out what fits with number of columns
    def nextpow2(i):
        n = 1
        while 2**n < i: n += 1 
        return n
    Nfft = int(2**nextpow2(Leng))
    print('The length of image column is {} and NFFT is {}'.format(Leng,Nfft))
    fdiv = np.linspace(0,1,int((Nfft/2)))
    f = Fs*fdiv
    Mtf = fft(Lsf,Nfft)
    ynorm = 2.0/N*np.abs(Mtf[0:Nfft//2])
    yfMax = max(ynorm)
    i = 0
    while i < len(ynorm):
        ynorm[i] = ynorm[i]/yfMax
        i = i+1    

    #now smooth the y data
    #yforplot = smooth(ynorm,19)
    yforplot = signal.savgol_filter(ynorm,19,3)

    print('Lengths of f = {} | Mtf = {} | and normalized MTF = {}'.format(len(f),len(Mtf),len(yforplot)))

    return f, xvals, dxvals,normESF,Lsf,ynorm,yforplot,ImgROI

currentPath = (os.path.dirname(os.path.realpath(__file__)))
InputPath = currentPath+'/ExampleImagesForMTF/Thermal450um_2-14-22/'
InputName = '10s_3em7_rotated.tif'
Input = InputPath+InputName
print('Input path is {}'.format(Input))
EffectivePixelSize = 0.102  #mm/pixel  for all except when using focal length extender
#EffectivePixelSize = 0.068  #mm/pixel this is the effective pixel size when using the focal length extender
MTF_output = CalculateMTF(Input,EffectivePixelSize,[0,0,0,0],RoiNeeded=True)

InputPath = currentPath+'/ImagesForMTF/Fast6-24-22/Filtered/'
InputName = 'Edge6mmTa_30s30EM_6.tif'
Input = InputPath+InputName
MTF_output2 = CalculateMTF(Input,EffectivePixelSize,[0,0,0,0],RoiNeeded=True)
#now make a 10% MTF dotted line
x10 = np.ones(len(MTF_output[0]))
x10 = x10*0.1
""" uncomment this section if you want to save the ESF,LSF, and MTF as .csv
ESF_array = np.array([MTF_output2[1],MTF_output2[3]])
np.savetxt(InputPath+'ESF_outputs.csv',np.transpose(ESF_array),delimiter=',')

LSF_array = np.array([MTF_output2[2],MTF_output2[4]])
np.savetxt(InputPath+'LSF_outputs.csv',np.transpose(LSF_array),delimiter=',')

MTF_array = np.array([MTF_output2[0],MTF_output2[5],MTF_output[6]])
np.savetxt(InputPath+'MTF_outputs.csv',np.transpose(MTF_array),delimiter=',')
"""

plt.figure(1,dpi=120)
plt.plot(MTF_output[1],MTF_output[3])
plt.xlabel('Distance (Pixel)')
plt.ylabel('Intensity (Gray Value)')
plt.title('Edge Spread Function for 2mm Cd using LiF:Zns')
plt.figure(2,dpi=120)
plt.plot(MTF_output[2],MTF_output[4])
plt.xlabel('Distance (Pixel)')
plt.ylabel('Derivative of Intensity (dGv/dx)')
plt.figure(3,dpi=120)
plt.plot(MTF_output[0],MTF_output[6], label='MTF Smoothed - LiF:ZnS')
plt.plot(MTF_output2[0],MTF_output2[6], label='MTF Smoothed - PVT')
plt.plot(MTF_output[0],x10,label='10% MTF')
plt.legend()
plt.xlabel('Spatial Frequency (lp/mm)')
plt.ylabel('Modulation Transfer Factor')
plt.title('MTF for LiF:ZnS using Thermal Neutrons and PVT Scintillator Using Fast Neutrons')
plt.show()

#this closes the cv2 ROI selector once we close out of our plots
cv2.waitKey(1)
cv2.destroyAllWindows()