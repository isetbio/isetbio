function [coneIndices, conePositionsDegs, coneTypes, coneIndicesInSerializedList] = indicesForConesAtPositions(obj, targetPosDegs)
% Return the indices for cones at a list of positions
%
% Syntax:
%   [coneIndices, conePositionsDegs, coneTypes] = indicesForConesAtPositions(obj, targetPosDegs);
%
% Description:
%    Return the indices for all cones lying at the target positions 
%       (specified in degrees of visual angle)  along  with the actual
%       cone positions.
%
% Inputs:
%    obj           - The cone mosaic hex object
%    targetPosDegs - A matrix of [N x 2] of N target positions
%
% Outputs:
%    coneIndices       - an [N x 1] vector of indices of cones whose location is 
%       closest to the target positions.
%    conePositionsDegs - an [N x 2] vector of positions (in  visual degrees)
%       of  the corresponding cones.
%    coneTypes - an [N x 1] vector of the corresponding cone types
%
% Optional key/value pairs:
%    None.
%
% Example  usage:
%   
%   c = coneMosaicHex(7, 'fovDegs', 0.05)
%   c.visualizeGrid('ticksInVisualDegs', true);
%   targetPosDegs = [0 0; 0.01 0.01; 0.02 0.02];
%   [coneIndices, conePositions, coneTypes] = c.indicesForConesAtPositions(targetPosDegs);
%   theOI = ...
%   coneExcitations = c.compute(theOI);
%   targetConeExcitations = coneExcitations(coneIndices);
%
% History:
%    4/9/2019  NPC  ISETBIO TEAM, 2018

    idx = find(obj.pattern > 1);
    [iRows, iCols] = ind2sub(size(obj.pattern), idx); 
    
    xSupport = squeeze(obj.patternSupport(1, :, 1));
    ySupport = squeeze(obj.patternSupport(:, 1, 2));
    coneXcoords = xSupport(iCols);
    coneYcoords = ySupport(iRows);
    
    allConePositionsDegs = [coneXcoords(:) coneYcoords(:)]*1e6/obj.micronsPerDegree;
    [distancesFromTargets, coneIndicesInSerializedList] = pdist2(allConePositionsDegs, targetPosDegs, 'euclidean', 'Smallest', 1);
    coneIndicesInSerializedList = unique(coneIndicesInSerializedList);
    
    coneIndices = idx(coneIndicesInSerializedList);
    
    [iRows, iCols] = ind2sub(size(obj.pattern), coneIndices); 
    coneXcoordsDegs = xSupport(iCols)*1e6/obj.micronsPerDegree;
    coneYcoordsDegs = ySupport(iRows)*1e6/obj.micronsPerDegree;
    conePositionsDegs = [coneXcoordsDegs(:) coneYcoordsDegs(:)];
    coneTypes = obj.pattern(coneIndices);
    
    beVerbose = false;
    if (beVerbose)
        for k = 1:size(targetPosDegs,1)
            switch coneTypes(k)
                case 2
                   coneType = 'L-cone'; 
                case 3
                    coneType = 'M-cone'; 
                case 4
                    coneType = 'S-cone';
                otherwise
                    error('How can this be?')
            end

            fprintf('[%2.0f]: (%s) at (%2.3f %2.3f) is closest to target pos (%2.3f,%2.3f), distance: %2.4f degs\n', ...
                k, coneType,  conePositionsDegs(k,1), conePositionsDegs(k,2), ...
                targetPosDegs(k,1), targetPosDegs(k,2), distancesFromTargets(k));
        end
    end
    
    
end

