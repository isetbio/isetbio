function m = macularCreate(macDensity,wave)
% Returns a structure containing several measures of the macular pigment
%
%     m = macularCreate(macDensity,wave)
%
% The human retina contains a pigment that covers the central (macular)
% region. This macular pigment passes certain wavelengths of light more
% than others.  The pigment varies in density from central vision, where it
% is highest, to increasingly peripheral vision.
%  
% This function returns several measures of the macular pigment wavelength
% properties as a function of macular pigment density (high in the fovea,
% lower in the near fovea).
%
% The returned structure, t, includes a variety of derived terms. This
% should help to keep the relationship between entities straight.
%
% macDensity is the estimated (average) peak density of the pigment across
% a variety of observers.  They estimate the average (across observers)
% peak density to be 0.28, with a range of 0.17 to 0.48.
%
% m.name    'Convenient name'
% m.type    'macular'
% m.wave    
% m.unitDensity:   The spectral density function with a maximum value of 1.0
% m.density:       The density for this instance
%
% Useful formulae
%
%   Absorbance spectra are normalized to a peak value of 1.
%   Absorptance spectra are the proportion of quanta actually absorbed.
%   Equation: absorptanceSpectra = 1 - 10.^(-OD * absorbanceSpectra)
%
% The original macular densities values were taken from the Stockman site.
% Go to http://cvision.ucsd.edu, then click on Prereceptoral filters.  At
% this point in time, I think the Psychtoolbox and the new Stockman site
% are authoritative.
%
% The densities were derived by Sharpe and Stockman based on some data from
% Bone. The paper describing why they like these is in Vision Research; I
% have downloaded the paper to Vision Science/Reference PDF/cone
% sensitivities
%
% Examples:
%   m = macularCreate;
%
% Copyright ImagEval Consultants, LLC, 2013

%% 
if notDefined('macDensity'), macDensity = 0.35; end
if notDefined('wave'), wave = (400:700)'; end

m.name = 'default human macular pigment';
m.type = 'macular';
m.wave = wave;

% Read in the Sharpe macular pigment curve and normalize to unit density
% Typical peak macular density, Estimated by Sharpe in VR paper, 1999 is
% 0.28.  Yet, the data they provide are at 0.3521.  It is probably not
% important to return the unit density, but we do.
density  = ieReadSpectra('macularPigment.mat',wave);
m.unitDensity = density / 0.3521;

m.density = macDensity;

return


