function temporalEquivEccentricity = retinalLocationToTEE(locationThetaDegrees, locationRadius, leftOrRightEyeIn)
% 
% retinalLocationToTEE converts the coordinates of a location on the retina
% given in polar coordinates (theta, r) to temporalEquivEccentricity (TEE).
% The TEE is used to calculate the diameter of the receptive field size at
% that location, given by Fig. 5 Chichilnisky & Kalmar (2002). It
% compensates for the observation that "contour lines of constant RGC
% density are circular in the temporal half but elliptical in the nasal
% half of the retina" (Chichilnisky & Kalmar, pg. 2738, 2002).
% 
% TEE = sqrt((0.61*X^2)+Y^2) (corrected from the text)
% 
% Inputs:
%   locationTheta - the theta value of the location of the patch in polar 
%       coordinates, in um, normally from 0-15.
%   locationRadius - the radius value of the location of the patch in polar 
%       coordinates, in degrees, from 0-360.
%   leftOrRightEye - the eye is used to determine which angular values are
%       nasal or temporal; can be strings 'left' or 'right' or integers '0'
%       (left) or '1' (right).
% 
% Outputs:
%   temporalEquivEccentricity - the eccentricity (um) if the location in the
%       nasal half were mapped onto the temporal half according to lines of
%       constant RGC density. This is used to predict RGC density and RF size.
% 
% Example: 
% temporalEquivEcc = retinalLocationToTEE(90, 3, 'right');
% 
% 
% 9/2015 JRG 
% (c) isetbio



if strcmp(leftOrRightEyeIn,'right')
    leftOrRightEye = 0;
elseif strcmp(leftOrRightEyeIn,'left')
    leftOrRightEye = 1;
    %     elseif ~((leftOrRightEyeIn == 0) || (leftOrRightEyeIn == 1))
    %         leftOrRightEye = leftOrRightEyeIn;
    %     else
    %         error('Third input to retinalLocationToTEE must be 0, 1, or ''left''/''right''');
end


locationTheta = (pi/180)*locationThetaDegrees;

if ( (locationTheta > (pi/2) && locationTheta < (3*pi/2)) && leftOrRightEye==1 ) || ...
    ( (locationTheta < (pi/2) || locationTheta > (3*pi/2)) && leftOrRightEye==0 ) 
    
    [xrad, yrad] = pol2cart(locationTheta, locationRadius);
    
    aspectRatio = 0.61;
    temporalEquivEccentricity = sqrt((xrad*aspectRatio).^2 + yrad.^2);

else 
    
    temporalEquivEccentricity = locationRadius;
    
end
