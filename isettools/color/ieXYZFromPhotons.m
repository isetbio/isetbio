function XYZ = ieXYZFromPhotons(photons,wave)
%Convert photon spectral power distribution into CIE XYZ
%
%   XYZ = ieXYZFromPhotons(photons,wave)
%
%  Converts a spectral power distribution in photons to CIE XYZ values.
%  The routine converts photons into photons and then calls ieXYZFromphotons.
%  See the comments about units in that that routine.
%
%  The format for photons can be XW or RGB
%
%  See also:
%    ieXYZFromEnergy
%
% Copyright ImagEval Consultants, LLC, 2003.

% Programming Notes:
%
% A simple way to do this is to call Quanta2photons and ieXYZFromphotons,
% i.e.
%   XYZ = ieXYZFromEnergy(Quanta2Energy(wave,photons),wave);
% However, for large input size (e.g. 1k x 1k x 100). This method could be
% very slow. Another option is to compute the XYZ responsivity curve for
% Quanta units and apply that to input photons directly
%
%(HJ) ISETBIO TEAM, 2014

% Force data into XW format.
iFormat = vcGetImageFormat(photons, wave);
if isempty(iFormat), error('unknown input format'); end

switch iFormat
    case 'RGB'
        % RGB format
        [xwData,r,c] = RGB2XWFormat(photons);
    otherwise
        % XW format
        xwData = photons;
end

if size(xwData,2) ~= length(wave)
    error('Problem converting input variable photons into XW format.');
end

% The spectra of the photons are in the rows of xwData.  We read the XYZ
% color matching functions into the columns of S. The original data is
% measured for energy data
S = ieReadSpectra('XYZ', wave);

% Adjust S to work for quanta units
S = Quanta2Energy(wave, S')';

% Compute return value XYZ as three column matrix
if numel(wave) > 1,  dWave = wave(2) - wave(1);
else                 dWave = 10;   disp('10 nm band assumed');
end

XYZ = 683 * dWave * (xwData*S);

% If it was sent in RGB, return it in RGB
switch iFormat
    case 'RGB'
        XYZ = XW2RGBFormat(XYZ,r,c);
    otherwise
        % XW format, do nothing
end

end

