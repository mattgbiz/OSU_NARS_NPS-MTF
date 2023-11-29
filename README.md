# OSU_NARS_NPS-MTF
Scripts to calculated noise power spectrum (NPS) and modulation transfer function (MTF) from digital radiographs. These measures are useful in determining image acquisition quality.

ReadMe for Ohio State University Nuclear Analysis and Radiation Sensor (NARS) Codes for Noise Power Spectrum (NPS) and Modulation Transfer Function (MTF)
*******************************************************************************************************
Original Author: Matthew Bisbee
Affiliations: Ohio State University Dept. of Mechanical and Aerospace Engineering, Nuclear Engineering
	      Nuclear Analysis and Radiation Sensor (NARS) Laboratory
	      DOE NEUP Fellowship FY19
	      Points of Contact: Advisor Dr. Raymond Cao - cao.152@osu.edu
			         Author Matthew Bisbee - bisbee.11@osu.edu
********************************************************************************************************

Python script contents include: CalculateNPS.py, CreateNoiseMap.py, and MTFCalculationfromImage.py

Two folders contain sample data that was used in the calculation of noise power spectrum (NPS) and modluation transfer function (MTF) for the OSU fast and thermal neutron imaging instrument. 

General Information: 

Each python script is intended to run in Python3.X. I used Python3.8 and Python3.9 on a Windows 10 device. The scripts attempt to find images based from the current path to the script itself (note this may only work for Windows devices, finding the current path for a MAC or Linux device may be different). 

NPS and MTF are useful functions that are dependant on spatial frequency (f) typically in units of lp/mm. One can use NPS and MTF in calculation of detective quantum efficiency (DQE) DQE is a unitless parameter that is defined as signal to noise ratio (SNR) output divided by SNR input. If a system is an ideal detection system, DQE would be 1. The lowest possible DQE would be 0 which in which there is no output signal. To quantify DQE in a real detection system the following equation is created. DQE(f) = [(mu)*MTF(f)^2]/[NPS(f)*N] where mu is the average pixel intensity and N is the integrated intensity of the beam in units of photon/mm^2. The units of NPS are mm^2 so these cancel out to provide a dimensionless DQE. In this work, mu was taken as a normalized number between 0 and 1 as "gray value" of the image would be dependent on the detector 
(for example an 8-bit camera would have a maximum value of 255 while a 16-bit camera  would be 65535). For this reason I took it as transmitted intensity of an open beam so an average value just under 1 due to photon statistics. 

*******************************************************************************************************
NPS Script:

The NPS script is built to read the current working path of the .py script and then find the image relative to that path. The ImgPath and ImageName variables can be changed for new images but a few example images have been provided to show how it works. The user selects a rectangular region of interest (ROI) of the open beam image. The ROI is then subdivided into sub regions each 32x32 so the ROI needs to be at least 32 pixels in both directions (ideally it will be much larger than 32 more like at least 150 pixels of open beam). px and py are the effective pixel size (mm/pixel) for the imaging system so this will need to be updated based on the user's system.

The code goes through a while loop to calculate each sub region's mean, density, and fast Fourier Transform values then sum this up at each step. Finnally, NPS is calculated taking this frequency domain Sumfft and normalizing it to px,px,Nx,Ny,Mm. The output 2D map is reordered for more intellible viewing by having the low spatial frequency in the center of the image going to higher spatial frequency towards the edges. The one-dimensonal average of the two-dimensional image is calculated along the horizontal and vertical directions near but not including the very center pixel. The one-dimensional array is output as a .csv file for further plotting ubt both the 2D and 1D plots are provided for quick viewing. 

*******************************************************************************************************
MTF Script: 

The MTF scipt has the calculation of MTF put into a user-written function to make it easier to calculate and compare multiple images simultaneously. An example is provided comparing the OSU imaging steup with 1 cm thick PVT high light yield fast neutron scintillator and 450 um thick LiF:ZnS thermal neutron scintillator. The MTF function takes inputs of image path, effective pixel size, and region of interest coordinates (as well as a boolean option for if the ROI needs to be pulled from the image itself). Various pertinent parameters are calculated from the effective pixel size. The edge spread function (ESF) is the average pixel values going from open beam to attenuated edge region averaged along all rows of ROI. Line spread function is the derivative of the ESF with respect to the columns along the ESF. The LSF is input into the fast Fourier Transform function along with parameter Nfft. The output MTF is normalized to 1 and smoothed with a savgol filter to provide the plotted curve. The section of block comments can be uncommented to save the ESF,LSF, and MTF arrays as .csv otherwise, the code simply plots the MTF results without saving.
