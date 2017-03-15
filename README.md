# TFMatlab
Traction Force Microscopy analysis software for MATLAB

These MATLAB functions can be used to calculate Traction Force Maps from fluorescent images of tracker particles imbedded in a elastic substrate. The code is intended for use with experiments modeled after those presented in Butler et al. (2002) Trends Cell Biol. 12:79-84.

Currently, the software is deisgned to accept Nikon NIS elements formated ND2 image sequences. Images can be analyzed using either Particle Tracking Velocimitry or Particle Image Velocimetry

# Usage
Initial data processing is handled by either CalcualteTFM(), which uses PTV to determine particle displacement, or using CalculateTFM_piv(), which uses PIV.

## CalculateTFM***()
If the functions are called without any arguments, the user will be presented with a series of GUI dialogs in which they can specify what data to process and what parameters to use.

Arguments should be passed to the functions as Name, value pairs.

The resulting data is saved to an uncompressed *.mat file (-v7.3) and optionally returned by the function as a structure.

Type "help CalculateTFM" at the MATLAB prompt for details.

## TFMForceViewer()
Use this function to view previously processed TFM data. The resulting GUI tool can be used to step through individual frames of the sequence, and also save a movie.


# Citing TFMatlab
TFMatlab was initially presented in the suplementary materials of Kovari et al. (2016) Biophys J. 111:12 DOI:10.1016/j.bpj.2016.11.009. If you find this software helpful, please consider citing that paper.
