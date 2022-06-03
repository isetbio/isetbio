function nearbyRGCindices = neihboringRGCindices(obj, theRGCindex)
       
    maxNeighborsNum = obj.wiringParams.maxNeighborsNum;
    maxNormDistance = obj.wiringParams.maxNeighborNormDistance;

    % Find the indices of its neirboring RGCs
    [distancesToNearbyRGCs, nearbyRGCindices] = RGCconnector.pdist2(...
            obj.RGCRFcentroidsFromInputs, ...
            obj.RGCRFcentroidsFromInputs(theRGCindex,:), '', ...
            'smallest', maxNeighborsNum+1);
                    
    % Exclude the theRGCindex, which has 0 distance
    [nearbyRGCindices,ia] = setdiff(nearbyRGCindices, theRGCindex);
    distancesToNearbyRGCs = distancesToNearbyRGCs(ia);

    % Exclude nearbyRGCs that are further than a maxDistance
    if (isinf(obj.localRGCRFspacingsMicrons(theRGCindex)))
        error('local spacing should not be inf here')
    end
    maxDistance = maxNormDistance * obj.localRGCRFspacingsMicrons(theRGCindex);
    
    idx = find(distancesToNearbyRGCs <= maxDistance);
    nearbyRGCindices = nearbyRGCindices(idx);
    distancesToNearbyRGCs = distancesToNearbyRGCs(idx);
    
    % Sort them according to distance to the RGCindex
    [~,idx] = sort(distancesToNearbyRGCs, 'ascend');
    nearbyRGCindices = nearbyRGCindices(idx);
end