function microSaccadeTargetLikelihoodSpatialMap = ...
    computeMicrosaccadeTargetLikelihoodSpatialMap(obj, ...
    heatMapXYsupportArcMin)
% Method to compute the microsaccade target likelihood map (a 2D Gaussian)
%
% Syntax:
%   microSaccadeTargetLikelihoodSpatialMap = ...
%       computeMicrosaccadeTargetLikelihoodSpatialMap(obj, ...
%       heatMapXYsupportArcMin)
%
% Description:
%    This method computes the microsaccade target likelihood map, which is a
%    two dimensional Gaussian. The map computed, 
%    obj.microSaccadeTargetLikelihoodSpatialMap, is combined with the 
%    constantly updated (by the recent history of the eye movement path), 
%    heatmap, obj.heatMapWeight, to find the next micro-saccade target.
%        
%
% Inputs:
%    obj                                    - Object. A fixationalEM object
%    heatMapXYsupportArcMin                 - Vector. The XY Support in
%                                             arcmin of the heat map.
%
% Outputs:
%    microSaccadeTargetLikelihoodSpatialMap - Matrix. The generated target
%                                             likelihood spatial map for
%                                             the microsaccades.
%
% Optional key/value pairs:
%    None.
%
% History:
%    01/03/18  NPC  ISETBIO Team, 2018
%    05/15/18  jnm  Formatting
%    05/24/18  NPC  Comments
%

% Gaussian with sigma obj.microaccadeTargetSigma
[x, y] = meshgrid(heatMapXYsupportArcMin, heatMapXYsupportArcMin);
microSaccadeTargetLikelihoodSpatialMap = exp(-0.5 * ...
    (x / obj.fixationMapSpaceConstantArcMin) .^ 2) .* ...
    exp(-0.5 * (y / obj.fixationMapSpaceConstantArcMin) .^ 2);
microSaccadeTargetLikelihoodSpatialMap = ...
    microSaccadeTargetLikelihoodSpatialMap / ...
    max(microSaccadeTargetLikelihoodSpatialMap(:));

end