% UTILITY
%
% These perform various little useful computations. Some are passed the wvf
% structure, but none set its fields directly or access its fields other
% than via wvfGet.
%
% Files
%   wvfComputeConePSF              - Return cone PSF and cone SCE Fraction using a wavefront object with PSF
%   wvfDefocusDioptersToMicrons    - Convert defocus in diopters to defocus in microns
%   wvfDefocusMicronsToDiopters    - Convert defocus in wavefront microns to diopters.
%   wvfKeySynonyms                 - Convert cell array of key value pairs to canonical form
%   wvfLCAFromWavelengthDifference - Longitudinal chromatic aberration (LCA), in diopters, between wavelengths
%   wvfLoadThibosVirtualEyes       - Load Thibos mean and covariance of Zernicke coefficients for human data
%   wvfOSAIndexToVectorIndex       - Convert a list of OSA j values to a WVF toolbox index
%   wvfOSAIndexToZernikeNM         - Convert OSA single Zernike index to the Zernike 2-index standard indexing
%   wvfWave2idx                    - Convert the wavelength list (wList) to indices relative to the list
%   wvfZernikeNMToOSAIndex         - Convert from Zernike 2-index standard index to OSA single-zernike index
