function [apertureCoverage, geometricCoverage] =  retinalCoverage(obj)
% Compute aperture and geometric retinal co verage
%
% Syntax:
%   [apertureCoverage, geometricCoverage] = obj.retinalCoverage()
%
% Description:
%    This function computes the aperture and geometric retinal coverages.
%
% Inputs:
%    obj               - The cone mosaic object
%
% Outputs:
%    apertureCoverage  - The aperture coverage for the cone mosaic object
%    geometricCoverage - The geometric coverage for the cone mosaic object
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/17  NPC  ISETBIO Team, 2017
%    02/19/18  jnm  Formatting

    if isa(obj, 'coneMosaicHex')
        minXY = min(obj.coneLocsHexGrid, [], 1);
        maxXY = max(obj.coneLocsHexGrid, [], 1);
        conesNum = size(obj.coneLocsHexGrid, 1);
    else
        minXY = min(obj.coneLocs, [], 1);
        maxXY = max(obj.coneLocs, [], 1);
        conesNum = size(obj.coneLocs, 1);
    end

    % Use the cone locs to compute the area of retina that the mosaic
    % extends over.
    retinalWidth = (maxXY(1) - minXY(1) + obj.pigment.width);
    retinalHeight = (maxXY(2) - minXY(2) + obj.pigment.width);
    retinalArea = retinalWidth * retinalHeight;

    % Compute cone aperture area from pigment properties
    coneApertureArea = obj.pigment.pdWidth * obj.pigment.pdHeight;

    % Compute geometric cone area from pigment properties
    coneArea = obj.pigment.width * obj.pigment.height;

    % compute aperture coverage

    if (obj.eccBasedConeQuantalEfficiency)
        apertureDiametersMeters = obj.computeApertureDiameters() * 1e-6;
        apertureCollectingAreas = pi*(apertureDiametersMeters/2).^2;
        retinalAreaCoveredWithApertures = sum(apertureCollectingAreas);
    else
        retinalAreaCoveredWithApertures = conesNum * coneApertureArea;
    end
    apertureCoverage = retinalAreaCoveredWithApertures / retinalArea;

    % compute geometric coverage
    geometricCoverage = apertureCoverage / coneApertureArea * coneArea;
end


