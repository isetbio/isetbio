function lut = ieLUTInvert(inLUT, nSteps)
%% Calculate inverse lookup table (lut) at certain sampling steps
%
%    lut = ieLUTInvert(inLUT,resolution)
%
% Inputs:
%   inLUT:      gamma table that converts linear DAC values to linear RGB.
%   nSteps:     samling steps, the returned gamma table is sampled at
%               the number of points specified by nSteps.  If it is not
%               passed in, the default is 2048.
%
% Outputs:
%   lut:  The returned lookup table
%   
% Example:
%   d = displayCreate('OLED-Sony');
%   inLUT = displayGet(d, 'gamma');
%   lut = ieLUTInvert(inLUT, 2);
%   vcNewGraphWin; plot(lut)
%
% See also:
%   ieLUTDigital, ieLUTLinear
%
% (c) Imageval Consulting, LLC 2013
%
% 1/7/15  dhb  Changed convention for passed resolution to be the number of
%              samples (nSteps) in the returned table.

%% Check inputs
if notDefined('inLUT'), error('input lut required'); end
if notDefined('nSteps'), nSteps = 2048; end

%% Computes inverse gamma table
%  Loop over primaries
nInSteps = size(inLUT,1);
y  = 1 : nInSteps;
iY = linspace(0,(nSteps-1)/nSteps,nSteps);
lut = zeros(length(iY), size(inLUT, 2));
for ii = 1 : size(inLUT, 2)
    % sort inLUT, theoretically, inLUT should be monochrome increasing, but
    % sometimes, the intensity at very low light levels cannot be measured
    % and we just set all of them to 0
    [x, indx] = unique(inLUT(:, ii));
    lut(:, ii) = interp1(x, y(indx), iY(:), 'pchip');
    lut = ieClip(lut, 0, nSteps - 1);
end

end