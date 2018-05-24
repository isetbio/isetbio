function targetPosArcMin = ...
    posOfMaxLikelihood(obj, targetLikelihoodSpatialMap)
% Calculate the position of max likelihood
%
% Syntax:
%   targetPosArcMin = posOfMaxLikelihood(obj, targetLikelihoodSpatialMap)
%   targetPosArcMin = obj.posOfMaxLikelihood(targetLikelihoodSpatialMap)
%
% Description:
%    Calculate the position of maximum likelihood using the provided data.
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

[~, idx] = max(targetLikelihoodSpatialMap(:));
if (numel(idx) > 1)
    % More than one locations have the same max. Pick one at random
    index2 = randperm(numel(idx));
    idx = idx(index2(1));
end
[row, col] = ind2sub(size(targetLikelihoodSpatialMap), idx);
targetPosArcMin = [obj.heatMapXYsupport(col) obj.heatMapXYsupport(row)];

end