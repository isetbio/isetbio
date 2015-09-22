function plotPatchEccentricity(retinalTheta, retinalRadius, leftOrRightEye, temporalEquivEcc)
% plotPatchEccentricity: a utility function of the rgc class that creates a
% contour plot showing the temporal equivalent eccentricity overlaid on the
% retina, with a red 'x' marking the location of the patch.
% 
% Inputs:
% 
% Outputs:
% 
% Example:
% 
% 09/2015 JRG

% Build contour plot of TEE using retinalLocationtoTEE
theta = 0:0.1:2*pi+0.1;
rad = 0:0.2:15;
theta_gr = repmat(theta,1,length(rad));
rad_gr = repmat(rad,1,length(theta)); 
for i = 1:length(theta_gr)
    tee(i) = retinalLocationToTEE((180/pi)*theta_gr(i), rad_gr(i), leftOrRightEye);
end
[xrad, yrad] = pol2cart(theta_gr, rad_gr);
xradrs = reshape(xrad, length(theta), length(rad));
yradrs = reshape(yrad, length(theta), length(rad));
tee_rs = reshape(tee, length(theta), length(rad));
figure; contourf(xradrs, yradrs, tee_rs, 16);
hold on;
h = polar(retinalTheta*(pi/180), retinalRadius,'xr');
set(h,'markers',16);
axis square;
colorbar;

if strcmp(leftOrRightEye, 'left');
    title(sprintf('%s eye, rad = %2.2f, theta = %2.2f, TEE = %2.2f\nnasal----------------temporal', leftOrRightEye, retinalRadius, retinalTheta, temporalEquivEcc),'fontsize',16);
else
    title(sprintf('%s eye, rad = %2.2f, theta = %2.2f, TEE = %2.2f\ntemporal----------------nasal', leftOrRightEye, retinalRadius, retinalTheta, temporalEquivEcc),'fontsize',16);

end