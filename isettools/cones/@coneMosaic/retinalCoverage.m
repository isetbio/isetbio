function [apertureCoverage, geometricCoverage] =  retinalCoverage(obj)
%RETINALCOVERAGE  Compute aperture and geometric retinal coverage
%
%   [apertureCoverage, geometricCoverage] = obj.retinalCoverage()
%
% NPC ISETBIO Team, 2017
%
    if isa(obj, 'coneMosaicHex')
        minXY = min(obj.coneLocsHexGrid,[],1);
        maxXY = max(obj.coneLocsHexGrid,[],1);
        conesNum = size(obj.coneLocsHexGrid,1);
    else
        minXY = min(obj.coneLocs,[],1);
        maxXY = max(obj.coneLocs,[],1);
        conesNum = size(obj.coneLocs,1);
    end
    
    % Use the cone locs to compute the area of retina that the mosaic extends over
    retinalWidth  = (maxXY(1)-minXY(1) + obj.pigment.width);
    retinalHeight = (maxXY(2)-minXY(2) + obj.pigment.width);
    retinalArea   = retinalWidth * retinalHeight;
    
    % Compute cone aperture area from pigment properties
    coneApertureArea = obj.pigment.pdWidth * obj.pigment.pdHeight;
    
    % Compute geometric cone area from pigment properties
    coneArea = obj.pigment.width * obj.pigment.height;
   
    % compute aperture coverage
    retinalAreaCoveredWithApertures = conesNum * coneApertureArea;
    apertureCoverage = retinalAreaCoveredWithApertures / retinalArea;
    
    % compute geometric coverage
    retinalAreaCoveredWithCones = conesNum * coneArea;
    geometricCoverage = retinalAreaCoveredWithCones / retinalArea;
end

