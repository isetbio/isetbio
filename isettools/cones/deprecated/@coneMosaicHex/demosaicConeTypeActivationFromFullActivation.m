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
    p.addParameter('interpolationMethod', 'linear', @validateInterpolationMethod);
    p.parse(varargin{:});
    interpolationMethod = p.Results.interpolationMethod;

    xPatternSupport = squeeze(obj.patternSupport(:,:,1));
    yPatternSupport = squeeze(obj.patternSupport(:,:,2));
    
    switch coneType
        case 'L-cones'
            targetConesIndices = find(obj.pattern==2);
            if (isstruct(interpolationMethod))
                kernelSigmaMicrons = interpolationMethod.kernelSigmaMicrons(1);
            end
        case 'M-cones'
            targetConesIndices = find(obj.pattern==3);
            if (isstruct(interpolationMethod))
                kernelSigmaMicrons = interpolationMethod.kernelSigmaMicrons(2);
            end
        case 'S-cones'
            targetConesIndices = find(obj.pattern==4);
            if (isstruct(interpolationMethod))
                kernelSigmaMicrons = interpolationMethod.kernelSigmaMicrons(3);
            end
        case 'All-cones'
            targetConesIndices = find(obj.pattern > 1);
            if (isstruct(interpolationMethod))
                kernelSigmaMicrons = mean(interpolationMethod.kernelSigmaMicrons);
            end
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
        
    % Compute demosaic response using griddata interpolation methods
    if (ischar(interpolationMethod))
        demosaicedResponseMap = griddata(coneXlocsDegs, coneYlocsDegs, coneResponses(:), X,Y, interpolationMethod);
    else
        % Linear griddata
        demosaicedResponseMap = griddata(coneXlocsDegs, coneYlocsDegs, coneResponses(:), X,Y, 'linear');
        
        % Find any samples that are nan (near the edges) and fill them with
        % the mean over the non-nan samples
        idx = find(~isnan(demosaicedResponseMap));
        meanD = mean(demosaicedResponseMap(idx));
        idx = find(isnan(demosaicedResponseMap));
        demosaicedResponseMap(idx) = meanD;

        % Generate smoothing kernel
        kernelSigmaDegs = kernelSigmaMicrons/obj.micronsPerDegree;
        kernelSigmaSamples = kernelSigmaDegs/degIncrement;
        kernelSizeSamples = round(kernelSigmaSamples*3)*2+1;
        fprintf('%s, kernel width: %d (%d)\n', coneType, kernelSizeSamples, numel(spatialSupportDegs));
        kernel = fspecial('gaussian',kernelSizeSamples,kernelSigmaSamples);
        
        % Add some extra margin pixels and fill them with the mean value
        margin = round(kernelSigmaSamples*3)+1;
        extendedSize = numel(spatialSupportDegs)+margin;
        originalSize = size(demosaicedResponseMap);
        extendedDemosaicedResponseMap = zeros(extendedSize,extendedSize) + meanD;
        
        % Paste the original values in the central part
        extendedDemosaicedResponseMap(margin+(1:originalSize(2)), margin+(1:originalSize(1))) = demosaicedResponseMap;
        % Convolve with kernel
        extendedDemosaicedResponseMap = conv2(extendedDemosaicedResponseMap, kernel, 'same');
        % back to original size
        demosaicedResponseMap = extendedDemosaicedResponseMap(margin+(1:originalSize(2)), margin+(1:originalSize(1)));
    end
    
end

function vStatus = validateInterpolationMethod(x)
    vStatus = true;
    
    % Check 1. Input is either a struct or a string
    if (~isstruct(x)) && (~ischar(x))
        error('Passed interolation method is neither a struct nor a string.');
    end
    
    % Check 2. Input is char. Valid values are 'linear' or 'nearest'
    if (ischar(x))
        if (~strcmp(x, 'linear')) && (~strcmp(x, 'nearest'))
            error('Passed interolation method is neither ''linear'' nor ''nearest''.');
        end
    end
    
    % If we made it here, the interpolation method is a struct
end
