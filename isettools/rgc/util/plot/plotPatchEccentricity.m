function plotPatchEccentricity(retinalTheta, retinalRadius, ...
    leftOrRightEye, temporalEquivEcc)
% Create a contour plot for temporal equivalent eccentricity on the retina
%
% Syntax:
%   plotPatchEccentricity(retinalTheta, retinalRadius, ...
%       leftOrRightEye, temporalEquivEcc)
%
% Description:
%    A utility function of the rgc class that creates a contour plot
%    showing the temporal equivalent eccentricity overlaid on the retina,
%    with a red 'x' marking the location of the patch.
%
%    This function contains examples of usage inline. To access these, type
%    'edit plotPatchEccentricity.m' into the Command Window.
%
% Inputs:
%    retinalTheta     - Numeric. The retinal angle.
%    retinalRadius    - Numeric. The retinal radius.
%    leftOrRightEye   - String. The eye side. Options are right and left.
%    temporalEquivEcc - Numeric. The temporal equivalent eccentricity.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * See retinalLocationToTEE for the details.
%    * Fig. 5, Chichilnisky & Kalmar, pg. 2738, 2002).
%           TEE = sqrt((0.61*X^2)+Y^2) (corrected from the text)
%

% History:
%    09/XX/15  JRG  Created
%    06/06/19  JNM  Documentation pass

% Examples:
%{
    plotPatchEccentricity(30, 5, 'right', 1)
%}

% Build contour plot of TEE using retinalLocationtoTEE
% Set up r, theta independent variables
theta = 0:0.1:2 * pi + 0.1;
rad = 0:0.2:15;
theta_gr = repmat(theta, 1, length(rad));
rad_gr = repmat(rad, 1, length(theta)); 
% find TEE at every point of r, theta grid
for i = 1:length(theta_gr)
    tee(i) = retinalLocationToTEE((180 / pi) * theta_gr(i), ...
        rad_gr(i), leftOrRightEye);
end

% Plot TEE
[xrad, yrad] = pol2cart(theta_gr, rad_gr);
xradrs = reshape(xrad, length(theta), length(rad));
yradrs = reshape(yrad, length(theta), length(rad));
tee_rs = reshape(tee, length(theta), length(rad));

figure; 
hold on;
contourf(xradrs, yradrs, tee_rs, 16);
h = polar(retinalTheta * (pi / 180), retinalRadius, 'xr');
set(h, 'markers', 16);
axis square;

ch1 = colorbar;
ylabel(ch1, 'TEE (mm)', 'fontsize', 18);
xlabel('Distance from fovea (mm)');
ylabel('Distance from fovea (mm)');

if strcmp(leftOrRightEye, 'left')
    title(sprintf(strcat('%s eye, rad = %2.2f, theta = %2.2f, ', ...
        'TEE = %2.2f\nnasal----------------temporal'), leftOrRightEye, ...
        retinalRadius, retinalTheta, temporalEquivEcc));
else
    title(sprintf(strcat('%s eye, rad = %2.2f, theta = %2.2f, ', ...
        'TEE = %2.2f\ntemporal----------------nasal'), leftOrRightEye, ...
        retinalRadius, retinalTheta, temporalEquivEcc));
end
set(gca, 'fontsize', 16);