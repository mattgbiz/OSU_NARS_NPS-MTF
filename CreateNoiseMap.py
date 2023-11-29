#The goal of this script is to create an array and save it as a tiff that is white noise with poisson distribution around mean lamb and pixel number x by y
import numpy as np
import cv2
lamb = 1000     #expected value of the poisson distribution (given by Padgett and Cao)
Nx = 1024       #sub-region size: Columns
Ny = 1024       #sub-region size: Rows

PoissonArray = np.random.poisson(lamb,(Nx,Ny))

#now I want to savet this array as a tif
Path = 'C:/Users/mattg_000/Documents/Research/ImageQuality_SNR_MTF_NPS_DQE/PoissonNoiseArray_2.tif'
cv2.imwrite(Path,PoissonArray)
