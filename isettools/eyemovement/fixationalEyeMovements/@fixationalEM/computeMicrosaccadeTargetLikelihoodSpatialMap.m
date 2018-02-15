% Method to compute the microsaccade target likelihood map (a 2D Gaussian)
% 
function microSaccadeTargetLikelihoodSpatialMap = computeMicrosaccadeTargetLikelihoodSpatialMap(obj, heatMapXYsupportArcMin)
    % Gaussian with sigma obj.microaccadeTargetSigma
    [x,y] = meshgrid(heatMapXYsupportArcMin,heatMapXYsupportArcMin);
    microSaccadeTargetLikelihoodSpatialMap = exp(-0.5*(x/obj.fixationMapSpaceConstantArcMin).^2) .* exp(-0.5*(y/obj.fixationMapSpaceConstantArcMin).^2);
    microSaccadeTargetLikelihoodSpatialMap = microSaccadeTargetLikelihoodSpatialMap / max(microSaccadeTargetLikelihoodSpatialMap(:));
end