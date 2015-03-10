function lens = lensCreate(lensDensity,wave, varargin)
% Returns a human lens structure
%
%     lens = lensCreate(lensDensity,wave)
%
% This function returns several measures of the lens in a structure.
%
% The returned structure, t, includes a variety of derived terms. This
% should help to keep the relationship between entities straight.
%
% macDensity is the estimated (average) peak density of the pigment across
% a variety of observers.  They estimate the average (across observers)
% peak density to be 0.28, with a range of 0.17 to 0.48.
%
% lens.name    'Convenient name'
% lens.type    'lens'
% lens.wave    
% lens.unitDensity:   The spectral density function (max =  1.0)
% lens.density:       The density for this instance
%
% Useful formulae
%
%   Absorbance spectra are normalized to a peak value of 1.
%   Absorptance spectra are the proportion of quanta actually absorbed.
%   Equation: absorptanceSpectra = 1 - 10.^(-OD * absorbanceSpectra)
%
% The original lens densities values were taken from PTB and/or the
% Stockman site. Go to http://cvision.ucsd.edu, then click on Prereceptoral
% filters.  At this point in time, I think the Psychtoolbox and the new
% Stockman site are authoritative.
%
% Examples:
%   lens = lensCreate;
%   
%
% HJ/BW ISETBIO Team 2013.

%% 
if notDefined('lensDensity'), lensDensity = 1; end
if notDefined('wave'),        wave = (400:700)'; end

lens.name = 'default human lens';
lens.type = 'lens';
lens.wave = wave;

% We need to figure out what the true density is.  It appears that the
% stored file is the true density, not the unit density.  If that is the
% case then the lens density here just scales the true density rather than
% providing the absolute density.
density          = ieReadSpectra('lensDensity.mat',wave);
lens.unitDensity = density;

lens.density     = lensDensity;

return


