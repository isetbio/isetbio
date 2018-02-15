function targetPosArcMin = posOfMaxLikelihood(obj,targetLikelihoodSpatialMap)
    [~,idx] = max(targetLikelihoodSpatialMap(:));
    if (numel(idx) > 1)
        % More than one locations have the same max. Pick one at random
        index2 = randperm(numel(idx));
        idx = idx(index2(1));
    end
    [row,col] = ind2sub(size(targetLikelihoodSpatialMap), idx);
    targetPosArcMin = [obj.heatMapXYsupport(col) obj.heatMapXYsupport(row)];
end