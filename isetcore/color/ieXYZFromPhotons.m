function XYZ = ieXYZFromPhotons(photons, wave)
% Convert spectral power distribution in photon units into CIE XYZ
%
% Syntax:
%   XYZ = ieXYZFromPhotons(photons, wave)
%
% Description:
%    Computes CIE XYZ values from a spectral power distribution in photon
%    units.
%
%    The format for photons can be XW or RGB
%
%    This function contains examples of usage. To acces, type 'edit
%    ieXYZFromPhotons.m' into the Command Window.
%
% Inputs:
%    photons - Vector. XW or RGB formatted spectral power distribution.
%    wave    - Vector. The wavelength(s).
%
% Outputs:
%    XYZ     - Vector. A CIE XYZ spectral power distribution, in the same
%              format (XW or RGB) as the input.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%   ieXYZFromEnergy, Quanta2Energy, vcGetImageFormat
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC.
%    xx/xx/14  HJ   ISETBIO TEAM, 2014
%    10/27/17  jnm  Comments & formatting
%    11/16/17  jnm  Formatting
%    07/10/19  JNM  Formatting update.

% Examples:
%{
   wave = 400:5:700;
   d65Energy = ieReadSpectra('D65', wave);
   d65Photons = Energy2Quanta(wave(:), d65Energy(:));
   d65XYZ = ieXYZFromPhotons(d65Photons(:)', wave)
%}

% Force data into XW format.
iFormat = vcGetImageFormat(photons, wave);
if isempty(iFormat), error('unknown input format'); end

switch iFormat
    case 'RGB'
        % RGB format
        [xwData, r, c] = RGB2XWFormat(photons);
    otherwise
        % XW format
        xwData = photons;
end

if size(xwData, 2) ~= length(wave)
    error('Problem converting input variable photons into XW format.');
end

% The spectra of the photons are in the rows of xwData. We read the XYZ
% color matching functions into the columns of S. The original data is
% measured for energy data
S = ieReadSpectra('XYZ.mat', wave);

% Adjust S to work for quanta units
S = Quanta2Energy(wave, S')';

% Compute return value XYZ as three column matrix
if numel(wave) > 1
    dWave = wave(2) - wave(1);
else
    dWave = 10;
    disp('10 nm band assumed');
end

XYZ = 683 * dWave * (xwData * S);

% If it was sent in RGB, return it in RGB
switch iFormat
    case 'RGB'
        XYZ = XW2RGBFormat(XYZ, r, c);
    otherwise
        % XW format, do nothing
end

end
