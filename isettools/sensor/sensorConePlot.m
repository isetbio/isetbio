function [support, spread, delta, coneMosaicImage] = sensorConePlot(sensor,support,spread,delta, whiteBackground)
% Plot the sensor cone array in saturated colors
%
%   [support, spread, delta, [coneMosaicImage]] = ...
%      sensorConePlot(sensor,[support],[spread],[delta], [whiteBackground])
%
% This should only be used for sensors based on a human cone-mosaic. This
% plotting routine will not run correctly (or at all) for typical ISET CMOS
% sensor mosaics.
% 
% The parameters support, spread and delta define the appearance of the
% plotted cone mosaic.  Each cone is rendered as a small blurry gaussian
% with support and spread as passed in.  The spatial sampling is delta. The
% whiteBackground parameters is a boolean indicating whether to generate an
% image with white or black background If a fourth output argument is
% present, the function returns the generated RGB image instead of plotting
% it.
%
% See also:  humanConeMosaic, conePlot, sensorCreateConeMosaic,
%            sensorShowCFA, sensorCreate
%
% Examples:
%
% All filled up, old defaults
%   sensor = sensorCreate('human');
%   sensorConePlot(sensor)
%
% A much smaller sensor with space in between
%   support = [7 7]; spread = 3; delta = .1;
%   sensorConePlot(sensor,support,spread,delta)
%
% Some empty spots (black).  New imaging defaults.
%   support = [5 5]; spread = 2; delta = .2;
%   sensorConePlot(sensor,support,spread,delta)
%
% Copyright ImagEval, 2010

if notDefined('support'), support = []; end
if notDefined('spread'), spread = []; end
if notDefined('delta'), delta = []; end
if notDefined('whiteBackground'), whiteBackground = false; end

% Read cone mosaic parameters and call cone plotting routine to make it
% look reasonably nice.
xy       = sensorGet(sensor,'human cone locs');
coneType = sensorGet(sensor,'cone type');
if isempty(xy)
    % A regular block array with human cones (rather than random mosaic).
    % Even so, we show the whole sensor CFA, as below for conePlot.
    fullArray = 1;
    sensorShowCFA(sensor,fullArray);
else
    if (nargout < 4)
        [support, spread, delta] = conePlot(xy, coneType, support, spread, delta, whiteBackground);
    else
        [support, spread, delta, coneMosaicImage] = conePlot(xy, coneType, support, spread, delta, whiteBackground);
    end
end

end