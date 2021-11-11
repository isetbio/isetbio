function cropMosaicDataForDesiredROI(obj)
    % Crop data for desired ROI
    diff = abs(bsxfun(@minus, obj.coneRFpositionsDegs, obj.eccentricityDegs));
    idxInsideROI = find((diff(:,1) <= 0.5*obj.sizeDegs(1)) & (diff(:,2) <= 0.5*obj.sizeDegs(2)));
    
    % Update state
    obj.updateStateGivenKeptConeIndices(idxInsideROI);
end

