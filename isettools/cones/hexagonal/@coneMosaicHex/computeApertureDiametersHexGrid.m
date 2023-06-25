function apertureDiametersMicrons = computeApertureDiametersHexGrid(obj)
% Return the aperture diameters for serialized cones in an ecc-based mosaic
%
% Syntax:
%   apertureDiametersMicrons = computeApertureDiametersHexGrid(obj)
%
% Description:
%     Return the aperture diameters for serialized cones in an ecc-based mosaic
%     ordered so they correspond to the serialized 1D response as well as to
%     obj.coneLocsHexGrid and the obj.coneTypesHexGrid
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
    coneXYEccentricities = obj.coneLocsHexGrid;
    coneEccentricitiesInMeters = (sqrt(sum(coneXYEccentricities.^2,2)))';
    coneAnglesInDegrees = atan2(squeeze(coneXYEccentricities(:,2)), squeeze(coneXYEccentricities(:,1))) / pi * 180;
        
    [~, apertureMeters, ~] = coneSizeReadData(...
            'eccentricity',coneEccentricitiesInMeters,...
            'angle',coneAnglesInDegrees, ...
            'useParfor', obj.useParfor);
        
    apertureDiametersMicrons = diameterForCircularApertureFromWidthForSquareAperture(apertureMeters * 1e6);
end
