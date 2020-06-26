function TEE = retinalLocationToTEE(theta, rho, eyeSide)
% Temporal equivalent eccentricity for a visual field position
%
% Syntax:
%   TEE = retinalLocationToTEE(theta, rho, eyeSide)
%
% Description:
%    This function is used to build RGC density and RF size values at
%    different eccentricities.
%
%    The polar coordinates of a location on the retina are converted to
%    temporalEquivEccentricity (TEE).
%
%    This function compensates for the observation that "contour lines of
%    constant RGC density are circular in the temporal half but elliptical
%    in the nasal half of the retina"
%       (Chichilnisky & Kalmar, pg. 2738, 2002).
%
%    The TEE is used to calculate the diameter of the receptive field size
%    at that location, given by Fig. 5 Chichilnisky & Kalmar (2002).
%
%       TEE = sqrt((0.61 * X ^ 2) + Y ^ 2) (corrected from the publication)
%
%    This function contains examples of usage. To access these examples,
%    type 'edit retinalLocationToTEE.m' into the Command Window.
%
% Inputs:
%    rho     - Numeric. The distance from the fovea of the location of the
%              patch, in mm, normally from 5-15. 
%    theta   - Numeric. The angle (deg) of the location of the patch in polar
%              coordinates, where 0 is the x-axis 
%    eyeSide - String. The eye is used to determine which angular values
%              are nasal or temporal; can be strings 'left' or 'right' or
%              integers '0' (left) or '1' (right).
%
% Outputs:
%    TEE     - Numeric. The temporal Equivalent Eccentricity. The
%              eccentricity (mm) of the location in the nasal half is
%              mapped onto the equivalent temporal eccentricity. This is
%              equivalent with respect to RGC density.
%
% Notes:
%    * [Note: JNM - When we've moved beyond 2015 in support, change out
%       strcmp(lower(... for strcmpi(... in the eyeSide check calls.]
%

% History:
%    XX/XX/15  JRG/BW  ISETBIO Team, 2015
%    06/05/19  JNM     Documentation pass

% Example:
%{
    % Plot the equivalent eccentricity as a graph
    vcNewGraphWin;
    theta = 0:0.1:(2 * pi + 0.1);
    rho = 5 * ones(size(theta));
    [x, y] = pol2cart(theta, rho);
    plot(x, y, 'rx')
    ee = zeros(size(theta));

    for ii = 1:length(theta)
        d = ieRad2deg(theta(ii));
        ee(ii) = retinalLocationToTEE(d, 5, 'right');
    end
    x = x .* ee / 5;
    y = y .* ee / 5;
    hold on;
    plot(x, y, 'bo');
    axis equal;
    grid on
%}
%% Get eye side
if isnumeric(eyeSide)  % Do nothing, already numeric.
elseif strcmp(lower(eyeSide), 'right'), eyeSide = 0;
elseif strcmp(lower(eyeSide), 'left'), eyeSide = 1;
else, error('Incorrect eye side specification %s\n', eyeSide);
end

% Convert angle in degrees to radians
theta = (pi / 180) * theta;

%% Apply formula for TEE
if ((theta > (pi / 2) && theta < (3 * pi / 2)) && eyeSide == 1) || ...
        ((theta < (pi / 2) || theta > (3 * pi / 2)) && eyeSide == 0 )
    [xrad, yrad] = pol2cart(theta, rho);

    % Apparently some correction for other quadrants.  Here is the formula.
    aspectRatio = 0.61;
    TEE = sqrt((xrad * aspectRatio) .^ 2 + yrad .^ 2);
else  % Temporal side. No need for a correction.
    TEE = rho;
end

end