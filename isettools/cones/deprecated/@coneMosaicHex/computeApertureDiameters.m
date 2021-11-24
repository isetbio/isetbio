function apertureDiametersMicrons = computeApertureDiameters(obj)
% Return the aperture diameters for all cones in an ecc-based mosaic
%
% Syntax:
%   apertureDiametersMicrons = computeApertureDiameters(obj)
%
% Description:
%     Return the aperture diameters for all cones in an ecc-based mosaic
%
% Inputs:
%     The cone mosaic
%
% Outputs:
%    apertureDiametersMicrons  - aperture diameters for all cones
%
% Optional key/value pairs:
%    None.
%       
    coneIndices = find(obj.pattern > 1);
    coneXYEccentricities = obj.coneLocs(coneIndices,:) / obj.resamplingFactor;
    coneEccentricitiesInMeters = (sqrt(sum(coneXYEccentricities.^2,2)))';
    coneAnglesInDegrees = atan2(squeeze(coneXYEccentricities(:,2)), squeeze(coneXYEccentricities(:,1))) / pi * 180;
        
    [~, apertureMeters, ~] = coneSizeReadData(...
            'eccentricity',coneEccentricitiesInMeters,...
            'angle',coneAnglesInDegrees, ...
            'useParfor', obj.useParfor);
        
    apertureDiametersMicrons = diameterForCircularApertureFromWidthForSquareAperture(apertureMeters * 1e6);
end
