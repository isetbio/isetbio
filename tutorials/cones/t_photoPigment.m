%% t_photoPigment
%
% The photoPigment object represents the data needed to calculate the
% capture of light by the cone photopigment.
%
% To create this object with its default values use

pp = photoPigment

%% Wavelength
%
% The internal wavelength representation of the parameters is from 390:830
% at 1nm.  But the user selects their own representation of the wavelength
% manually.  By default this is 400:10:700.

%%  Absorbance
% The cone absorbance function used by default is read from the file
% data/human/coneAbsorbance
%
% It is stored in normalized form, by default,  at a 400:10:700 wavelength
% resolution
vcNewGraphWin; plot(pp.wave,pp.absorbance)


