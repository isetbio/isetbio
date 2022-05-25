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
    maxDistance = maxNormDistance * obj.RGCRFspacingsMicrons(theRGCindex);
    nearbyRGCindices = nearbyRGCindices(find(distancesToNearbyRGCs <= maxDistance));
end