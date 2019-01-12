function [demosaicedResponseMap, spatialSupportDegs, coneResponses, coneXlocsDegs, coneYlocsDegs] = ...
    demosaicConeTypeActivationFromFullActivation(obj, coneType,...
    theFullPatternResponse, demosaicingSampleSpacingMicrons, varargin)
%  Obtain a demosaiced map of the activation for a single cone type from
%  the full mosaic response. The spatial support for the demosaiced
%  response map is returned in degrees.  Also returns the non-demosaiced response 
%  for the specified cone type and the x/y coords of all cones of that cone type
%  
% Usage:
%
%  theFullMosaicIsomerizations = theConeMosaic.compute(theOI, 'currentFlag',false);
%  demosaicingSampleSpacingMicrons = 1.0;
%  coneTypeName = 'L-cones'; % choose from {'L-cones', 'M-cones', 'S-cones', 'All-cones'};
%  [demosaicedResponseMap, demosaicedMapSupportDegs, ...
%   theConeIsomerizations, xConeLocsDegs, yConeLocsDegs] = ...
%      theConeMosaic.demosaicConeTypeActivationFromFullActivation(coneTypeName, ...
%      theFullMosaicIsomerizations, demosaicingSampleSpacingMicrons);
%
%   figure()
%   imagesc(demosaicedMapSupportDegs, demosaicedMapSupportDegs,demosaicedResponseMap)
%
%
% NPC, ISETBIO TEAM, 2018

    % Parse optional input
    p = inputParser;
    p.addParameter('interpolationMethod', 'nearest', @ischar);
    p.parse(varargin{:});
    interpolationMethod = p.Results.interpolationMethod;

    xPatternSupport = squeeze(obj.patternSupport(:,:,1));
    yPatternSupport = squeeze(obj.patternSupport(:,:,2));
    
    switch coneType
        case 'L-cones'
            targetConesIndices = find(obj.pattern==2);
        case 'M-cones'
            targetConesIndices = find(obj.pattern==3);
        case 'S-cones'
            targetConesIndices = find(obj.pattern==4);
        case 'All-cones'
            targetConesIndices = find(obj.pattern > 1);
        otherwise
            error('Invalid cone type: ''%s''.', coneType);
    end
    
    if (isempty(targetConesIndices))
        demosaicedResponseMap = [];
        spatialSupportDegs = [];
        coneResponses = []; 
        coneXlocsDegs = [];
        coneYlocsDegs = [];
        return;
    end
    
    coneResponses = theFullPatternResponse(targetConesIndices);
    xPatternSupport = xPatternSupport(targetConesIndices);
    yPatternSupport = yPatternSupport(targetConesIndices);
    
    coneXlocsDegs = xPatternSupport*1e6/obj.micronsPerDegree;
    coneYlocsDegs = yPatternSupport*1e6/obj.micronsPerDegree;
        
    degIncrement = demosaicingSampleSpacingMicrons/obj.micronsPerDegree;
    
    demosaicRangeMicrons = max(obj.fov)*obj.micronsPerDegree;
    demosaicRangeDegs = demosaicRangeMicrons / obj.micronsPerDegree;
    spatialSupportDegs = (-demosaicRangeDegs/2):degIncrement:(demosaicRangeDegs/2);
    spatialSupportDegs = spatialSupportDegs - mean(spatialSupportDegs(:));
    [X,Y] = meshgrid(spatialSupportDegs,spatialSupportDegs);
        
    % Compute demosaic response
    demosaicedResponseMap = griddata(coneXlocsDegs, coneYlocsDegs, coneResponses(:), X,Y, interpolationMethod);
end