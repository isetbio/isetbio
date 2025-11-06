function visualizeOpticsAtTargetWavelength(theOI, targetWavelength, ...
    spatialSupportArcMin, spatialFrequencySupportCyclePerDeg, varargin)
% VISUALIZEOPTICSATTARGETWAVELENGTH Visualizes the optical properties at a specified target wavelength.
%
% Syntax:
%   visualizeOpticsAtTargetWavelength(theOI, targetWavelength, ...
%   spatialSupportArcMin, spatialFrequencySupportCyclePerDeg, varargin)
%
% Description:
%   This function visualizes the Point Spread Function (PSF) and Optical Transfer Function (OTF)
%   at a specified target wavelength. It creates a figure with two subplots: one for the PSF
%   and another for the OTF, allowing for analysis of optical performance at the given wavelength.
%
% Inputs:
%   theOI                     - The optical image object containing the optical properties.
%   targetWavelength          - The wavelength (in nanometers) at which to visualize the optics.
%   spatialSupportArcMin      - The spatial support in arc minutes for the PSF visualization.
%   spatialFrequencySupportCyclePerDeg - The spatial frequency support in cycles per degree for the OTF visualization.
%
% Optional Name-Value Pair Arguments:
%   'extraOTFData'           - Additional data for the OTF visualization (default: []), if needed.
%
% Example:
%   visualizeOpticsAtTargetWavelength(opticalImage, 550, 30, 10, 'extraOTFData', someData);
%
% See also: VISUALIZEPSF, VISUALIZEOTF

p = inputParser;
p.addParameter('extraOTFData', []);

% Parse input
p.parse(varargin{:});
extraOTFData = p.Results.extraOTFData;

hFig = figure(); clf;
set(hFig, 'Color', [1 1 1]);
ax1 = subplot(1,2,1);
visualizePSF(theOI, targetWavelength, spatialSupportArcMin, ...
    'axesHandle', ax1);

ax2 = subplot(1,2,2);
visualizeOTF(theOI, targetWavelength, spatialFrequencySupportCyclePerDeg, ...
    'axesHandle', ax2, 'extraData', extraOTFData);
end
