#the purpose of this code is to find the NPS for a given image
#the image should be flat fielded and cropped to the uniform beam area
#this means that it should be some square image of the open beam for the OSURR imaging station
#ideally one would crop the image to fit well with the sub-ROI chosen
#for example, an open beam that is 1024x1024 (like the poisson image) would be well suited...
    #for subregions that are 128x128 as that fits nicely within 1024x1024

import numpy as np
import cv2, os
from matplotlib import pyplot as plt

def split(array, nrows, ncols):
    """Split a matrix into sub-matrices."""
    r, h = array.shape
    return (array.reshape(h//nrows, nrows, -1, ncols)
                 .swapaxes(1, 2)
                 .reshape(-1, nrows, ncols))

currentPath = (os.path.dirname(os.path.realpath(__file__)))
ImgPath = currentPath+'/ExampleImagesForNPS/LiFZnS_450um10s_3em'
ImageName = '/White_image#1.tif'
OutputName = ImgPath+'/NPS_450umLiF_1.csv'
#read in the image
Img = cv2.imread(ImgPath+ImageName,-1)
ReadRealImage = True
if ReadRealImage == True:
    maximg = np.max(Img)
    #maximg = np.max(img_med1)
    imgmult = int(65530/maximg) #multiply this to stretch the image to fill 16bit depth for cv2
    cal = Img*imgmult
    roi = cv2.selectROI(cal)
    #print rectangle point of selected ROI
    roi_cropped = Img[int(roi[1]):int(roi[1]+roi[3]),int(roi[0]):int(roi[0]+roi[2])]
    imgWidth = roi[2]
    imgHeight = roi[3]
    #now eventually I may do some type of forcing width/height to be in some range like 256 or something
    #looks like the open beam image i have for around an effective pixel size of 0.1mm is 192x192 max which splits well into 32 subROI
    #the 192x192 was selected as larger than that the beam may not be uniform at the edges
    Nx = 32
    Ny = 32
    px = 0.102     #effective pixel size for specific radiograph
    py = 0.102     #effective pixel size for specific radiograph
    #px = 0.016      #actual dimensional size of pixel
    #py = 0.016      #actual dimensional size of pixel
else:
    #get the dimensions of the read in image
    imgSize = Img.shape
    imgWidth = int(imgSize[0])        
    imgHeight = int(imgSize[1])
    #as for the poisson tester, that was 1024x1024 and I'll use overlapping 128x128s
    #Nx = 128    #subregion width will be 32 for my images NEEDS TO BE EVEN
    #Ny = 128    #subregion height will be 32 for my images NEEDS TO BE EVEN
    #px = 0.2    #effective pixel width (mm) (will be whatever my effective pixel size is)
    #py = 0.2    #effective pixel height (mm)
    Nx = 32
    Ny = 32
    px = 0.1
    py = 0.1



Mm = 1      #this is the total number of regions summed which will be incremented each time
i = 0       #increments for the rows/columns to be considered
j = 0
#initialize a numpy array that is Ny by Nx in size
Sumfft = np.zeros([Ny,Nx])

while i < (imgHeight-(Ny/2)-1):
    while j < (imgWidth-(Nx/2)-1):
        #now for this region from i to i+Ny and from j to j+Nx we need the fft
        SubRegion = Img[i:i+Ny,j:j+Nx]
        #now with this subregion, find the mean value
        SubMean = np.mean(SubRegion)
        #for all values in the subregion, subtract off the mean (take the absolute value of these)
        Density = np.abs(SubRegion-SubMean)
        #now take the 2D fast fourier transform of the Density and add that to our total fft result
        #Partfft = np.square(np.real(np.fft.fft2(Density)))   #need to square this value as well: Also I am only taking the real component that is coming out from this COULD BE INCORRECT WAY OF DOING IT
        Partfft = np.square(np.real(np.log(np.fft.fft2(Density))))
        Sumfft = Sumfft+Partfft
        #increment j by whatever our subregion/2 is
        j = int(j+(Nx/2))
        Mm = Mm+1   #this is for how many times we went through for averaging afterwards
    #increment i by whatever our subregion/2 is and reset j
    i = int(i+(Ny/2))
    j = 0

#after leaving the loops, we have the final summed fft and the total number of iterations done
#now average this and normalize to the input parameters
Nps = Sumfft*px*py/Nx/Ny/Mm
#print(Nps)
#now I need to split the NPS into 4 quadrants and then rotate each of them by 180 degrees
#this is because the output has spatial frequency of 0 on the edges and I want those to be in the center
#then I can project this to 1 D
upperhalf = np.hsplit(np.vsplit(Nps,2)[0],2)
lowerhalf = np.hsplit(np.vsplit(Nps,2)[1],2)
Qtl = upperhalf[0]
Qtr = upperhalf[1]
Qbl = lowerhalf[0]
Qbr = lowerhalf[1]

#with these quadrants I can get the average horizontal/vertical arrays to average
horizontal1 = Qtl[0,:]
horizontal2 = np.flip(Qtr[0,:])
horizontal3 = Qbl[int((Ny/2)-1),:] #(Ny/2)-1,:
horizontal4 = np.flip(Qbr[int((Ny/2)-1),:])
vertical1 = Qtl[:,0]
vertical2 = Qtr[:,int((Ny/2)-1)]
vertical3 = np.flip(Qbl[:,0])
vertical4 = np.flip(Qbr[:,int((Ny/2)-1)])
pixels_per_mm = 1/px
Fs = np.linspace(0,pixels_per_mm,int(Ny/2))

plt.figure(1)
plt.plot(Fs[1:],horizontal1[1:],'k--',label='Horizontal TL')
plt.plot(Fs,horizontal2,'b--',label='Horizontal TR')
plt.plot(Fs,horizontal3,'r--',label='Horizontal BL')
plt.plot(Fs,horizontal4,'g--',label='Horizontal BR')
plt.plot(Fs[1:],vertical1[1:],'k-',label='Vertical TL')
plt.plot(Fs,vertical2,'b-',label='Vertical TR')
plt.plot(Fs,vertical3,'r-',label='Vertical BL')
plt.plot(Fs,vertical4,'g-',label='Vertical BR')
plt.legend()
plt.xlabel('Spatial Frequency (lp/mm)')
plt.ylabel('NPS (mm\^2)')
inc = 0
Average1D_NPS = np.zeros(len(horizontal1))
while inc < len(horizontal1):
    if inc == 0:
        Average1D_NPS[inc] = (horizontal2[inc]+horizontal3[inc]+horizontal4[inc]+vertical2[inc]+vertical3[inc]+vertical4[inc])/6
    else:
        Average1D_NPS[inc] = (horizontal1[inc]+horizontal2[inc]+horizontal3[inc]+horizontal4[inc]+vertical1[inc]+vertical2[inc]+vertical3[inc]+vertical4[inc])/8
    inc +=1

#save these as .csv so that I can plot with MTF data
NPS_array = np.array([horizontal1,horizontal2,horizontal3,horizontal4,vertical1,vertical1,vertical3,vertical4])
np.savetxt(ImgPath+'NPSlinesBeforeAvg_200umLiF_SwitchMedian_1.csv',np.transpose(NPS_array),delimiter=',')
np.savetxt(OutputName,Average1D_NPS,delimiter=',')


plt.figure(2)
plt.semilogy(Fs,Average1D_NPS)
plt.xlabel('Spatial Frequency (lp/mm)')
plt.ylabel('NPS (mm^2)')
plt.title('1D approximation of NPS')
plt.show()


TL_flip = np.rot90(Qtl,2)   #rotated by 90 degree twice aka 180 degree
TR_flip = np.rot90(Qtr,2)
BL_flip = np.rot90(Qbl,2)
BR_flip = np.rot90(Qbr,2)
#this combines them back together into one image to view
FlipCombine = np.vstack([np.hstack([TL_flip,TR_flip]),np.hstack([BL_flip,BR_flip])])

print('Maximum value of the array = {} | minimum value of the array = {} | average value of the array = {}'.format(np.max(Nps),np.min(Nps),np.mean(Nps)))
#now I want to see what the Nps looks like 
NPSplot = plt.imshow(FlipCombine)
plt.colorbar(NPSplot)
plt.show()