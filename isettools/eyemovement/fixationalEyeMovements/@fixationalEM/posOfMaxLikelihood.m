function targetPosArcMin = ...
    posOfMaxLikelihood(obj, targetLikelihoodSpatialMap)
% Return the (x,y,) position corresponding to the max likelihood map.
%
% Syntax:
%   targetPosArcMin = posOfMaxLikelihood(obj, targetLikelihoodSpatialMap)
%   targetPosArcMin = obj.posOfMaxLikelihood(targetLikelihoodSpatialMap)
%
% Description:
%    Return the (x,y,) position corresponding to the max likelihood map,
%    which is the a combination of a zero-centered Gaussian and the heat
%    map based on previous em locations. This position is going to be the
%    next micro-saccade target.
%
% Inputs:
%    obj                        - Object. A fixationalEM object.
%    targetLikelihoodSpatialMap - Matrix. A matrix containing the spatial
%                                 map of the target likelihood.
%
% Outputs:
%    targetPosArcMin            - Vector. A 1x2 Vector containing the
%                                 most likely position.
%
% Optional key/value pairs:
%    None.
%
% History:
%    01/03/18  NPC  ISETBIO Team, 2018
%    05/15/18  jnm  Formatting
%    05/24/18  NPC  Comments

[~, idx] = max(targetLikelihoodSpatialMap(:));
if (numel(idx) > 1)
    % More than one locations have the same max. Pick one at random
    index2 = randperm(numel(idx));
    idx = idx(index2(1));
end
[row, col] = ind2sub(size(targetLikelihoodSpatialMap), idx);
targetPosArcMin = [obj.heatMapXYsupport(col) obj.heatMapXYsupport(row)];

end