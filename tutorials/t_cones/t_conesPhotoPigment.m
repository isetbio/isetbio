%%t_conesPhotoPigment  Illustrate photoPigment object
%
% Description:
%     The photoPigment object represents the data needed to calculate the
%     capture of light by cones.  This tutorial illustrates how it works.

%% Initialize
ieInit;

%% Create photoPigment object with default values
pp = photoPigment;

%% Wavelength
%
% The internal wavelength representation of the parameters is from 390:830
% at 1 nm.  The user can set the wavelength representation used for interface
% with other routines. By default this is 400:10:700.

%%  Absorbance
%
% The cone absorbance function used by default is obtained
% via data routine getLogConeAbsorbance and converted 
% to straight (non-log) absorbance when it is stored in the object.
%
% Absorbance is sometimes called optical density.
%
% The peak absorbance is 1 by convention in how the values
% are tabulated. In the plot below they are close to 1 because
% the actual peak occurs at a wavelength between the default
% 10 nm sampling steps.
vcNewGraphWin; plot(pp.wave,pp.absorbance)

%% Geometry
ecc = 0.0; ang = 0;

% Spacing (um), aperture (um), density (cones/mm2)
[s,a,d] = coneSize(ecc,ang)
