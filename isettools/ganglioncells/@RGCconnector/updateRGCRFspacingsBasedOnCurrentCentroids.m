function updateRGCRFspacingsBasedOnCurrentCentroids(obj)
    % Compute RGCRF spacings from their current centroids 
    obj.RGCRFspacingsMicrons = RGCmodels.Watson.convert.positionsToSpacings(obj.RGCRFcentroidsFromInputs);
end