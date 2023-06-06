function generateSamplingGrids(obj, minSpatialSamplingDegs)

    centerConnectableConeTypes = obj.retinalRFmodelParams.centerConnectableConeTypes;

    % Find out the range of cones in the RF center
    allConesNumPooledByTheRFcenters = full(sum(obj.theRGCMosaic.rgcRFcenterConeConnectivityMatrix,1));
    conesNumPooledByTheRFcenters = unique(allConesNumPooledByTheRFcenters);
    fprintf('Cones/RF center for this mosaic: %d\n', conesNumPooledByTheRFcenters);

    % Generate grid of spatial positions x conesNumInRFcenter
    [multiFocalPositionGrid, multiFocalConesNumInRFcenterGrid] = ...
        meshgrid(1:size(obj.nominalSpatialSamplingGrid,1), conesNumPooledByTheRFcenters);
    
    % Allocate memory
    multiFocalPositionsNum = numel(multiFocalConesNumInRFcenterGrid);
    spatialPositionsVisitedForConesNumPooled = cell(1, max(conesNumPooledByTheRFcenters));

    obj.conesNumPooledByTheRFcenterGrid = [];
    obj.visualSTFSurroundToCenterRcRatioGrid = [];
    obj.visualSTFSurroundToCenterIntegratedSensitivityRatioGrid = [];
    obj.targetRGCindicesWithLconeMajorityCenter = [];
    obj.targetRGCindicesWithMconeMajorityCenter = [];

    for iMultifocalPositionIndex = 1:multiFocalPositionsNum

        % Cones num in RF center
        conesNumPooled = multiFocalConesNumInRFcenterGrid(iMultifocalPositionIndex);

        % Sampling position (within the  mosaic)
        eccPositionIndex = multiFocalPositionGrid(iMultifocalPositionIndex);

        % Start with the nominal position
        nominalPositionDegs = obj.nominalSpatialSamplingGrid(eccPositionIndex,:);

        % Get indices of center cones for the RGC that is closest to the nominalPositionDegs
        indicesOfRGCsWithThisManyCenterCones = find(allConesNumPooledByTheRFcenters == conesNumPooled);
        distances = sqrt(sum((bsxfun(@minus, obj.theRGCMosaic.rgcRFpositionsDegs(indicesOfRGCsWithThisManyCenterCones,:), nominalPositionDegs)).^2,2));
        [~,idx] = sort(distances, 'ascend');
        indicesOfRGCsWithThisManyCenterCones = indicesOfRGCsWithThisManyCenterCones(idx);

   
        % Examine the nearby RGCs to determine the one with the
        % highest spectral purity (smallest ratio of cone type1 to cone
        % type2)
        ratio = zeros(1,numel(indicesOfRGCsWithThisManyCenterCones));
        for i = 1:numel(indicesOfRGCsWithThisManyCenterCones)
            theTargetRGCindex = indicesOfRGCsWithThisManyCenterCones(i);
            theCenterConeTypeWeights = obj.centerConeTypeWeights(theTargetRGCindex);
            ratio(i) = theCenterConeTypeWeights(2)/theCenterConeTypeWeights(1);
        end

        [~,idx] = min(ratio);
        theTargetRGCindex = indicesOfRGCsWithThisManyCenterCones(idx);
        [~, ~, theMajorityConeTypeOfTheTargetRGC] = obj.centerConeTypeWeights(theTargetRGCindex);

        % Get nearby RGC index with different center majority cone type
        [theTargetRGCindexOfDifferentMajorityConeType, theMajorityConeTypeOfTheTargetRGCWithDifferentMajorityConeType] = ...
            nearbyRGCindexOfDifferentMajorityCenterConeType(obj, theTargetRGCindex, indicesOfRGCsWithThisManyCenterCones);

        % Update the nominalPositionDegs for this RTVFobj to reflect the actual position of theTargetRGCindex
        nominalPositionDegs = obj.theRGCMosaic.rgcRFpositionsDegs(theTargetRGCindex,:);

        % Check to see if we visited a nearby position
        alreadyVisitedPositions = spatialPositionsVisitedForConesNumPooled{conesNumPooled};
        if (~isempty(alreadyVisitedPositions))
            mindistance = min(sqrt(sum((bsxfun(@minus, alreadyVisitedPositions, [nominalPositionDegs(1) nominalPositionDegs(2)])).^2,2)));
            if (mindistance < minSpatialSamplingDegs)
                continue;
            end
        end
       
        % We havent visited, so add it
        alreadyVisitedPositions = cat(1, alreadyVisitedPositions, [nominalPositionDegs(1) nominalPositionDegs(2)]);
        spatialPositionsVisitedForConesNumPooled{conesNumPooled} = alreadyVisitedPositions;

        indicesOfConesPooledByTheRFcenter = find(obj.theRGCMosaic.rgcRFcenterConeConnectivityMatrix(:,theTargetRGCindex)> 0);
        typesOfConesPooledByTheRFcenter = obj.theRGCMosaic.inputConeMosaic.coneTypes(indicesOfConesPooledByTheRFcenter);

        % Assert that we have the correct number of center cones
        assert((numel(indicesOfConesPooledByTheRFcenter) == conesNumPooled), ...
            sprintf('indicesOfConesPooledByTheRFcenter should have %d entries but it has %d.', numel(indicesOfConesPooledByTheRFcenter), conesNumPooled));
       
        % Assert that the center cones are all connectable
        assert(all(ismember(typesOfConesPooledByTheRFcenter, centerConnectableConeTypes)), ...
            sprintf('indicesOfConesPooledByTheRFcenter are not all connectable'));

        % Rs/Rc ratio
        surroundToCenterRcRatio = RGCmodels.CronerKaplan.constants.surroundToCenterRcRatio;

        % Temporal-equivalent eccentricity based SCint sensitivity ratio
        temporalEquivalentEccDegs = obj.theRGCMosaic.temporalEquivalentEccentricityForEccentricity(nominalPositionDegs);
        radialTemporalEquivalentEccDegs = sqrt(sum(temporalEquivalentEccDegs.^2,2));
        scIntSensitivity = RGCmodels.CronerKaplan.constants.surroundToCenterIntegratedSensitivityRatioFromEccDegsForPcells(radialTemporalEquivalentEccDegs);
       
        % The various grids
        obj.conesNumPooledByTheRFcenterGrid = ...
            cat(1, obj.conesNumPooledByTheRFcenterGrid,conesNumPooled);

        obj.visualSTFSurroundToCenterRcRatioGrid  = ...
            cat(1, obj.visualSTFSurroundToCenterRcRatioGrid, surroundToCenterRcRatio);

        obj.visualSTFSurroundToCenterIntegratedSensitivityRatioGrid = ...
            cat(1, obj.visualSTFSurroundToCenterIntegratedSensitivityRatioGrid, scIntSensitivity);

         % The target RGC indices
         switch (theMajorityConeTypeOfTheTargetRGC)
             case cMosaic.LCONE_ID
                 obj.targetRGCindicesWithLconeMajorityCenter = cat(1, ...
                     obj.targetRGCindicesWithLconeMajorityCenter, theTargetRGCindex);
                 obj.targetRGCindicesWithMconeMajorityCenter = cat(1, ...
                     obj.targetRGCindicesWithMconeMajorityCenter, theTargetRGCindexOfDifferentMajorityConeType);

             case cMosaic.MCONE_ID
                 obj.targetRGCindicesWithMconeMajorityCenter = cat(1, ...
                     obj.targetRGCindicesWithMconeMajorityCenter, theTargetRGCindex);
                 obj.targetRGCindicesWithLconeMajorityCenter = cat(1, ...
                     obj.targetRGCindicesWithLconeMajorityCenter, theTargetRGCindexOfDifferentMajorityConeType);
             otherwise
                 % equal ratios
                 error('equal ratios')
         end
    end % for iMultifocalPositionIndex 
end


function [theNearbyTargetRGCindexOfDifferentMajorityConeType, ...
          theNearbyMajorityConeType] = nearbyRGCindexOfDifferentMajorityCenterConeType(...
            obj, theTargetRGCindex, indicesOfRGCsWithThisManyCenterCones)

    [~, ~, theMajorityConeType] = obj.centerConeTypeWeights(theTargetRGCindex);
    notFound = true;
    iRGC = 1;
    while (notFound) && (iRGC < numel(indicesOfRGCsWithThisManyCenterCones))
        iRGC = iRGC + 1;
        theNearbyTargetRGCindex = indicesOfRGCsWithThisManyCenterCones(iRGC);
        [theCenterConeTypeWeights, ~, theNearbyMajorityConeType] = obj.centerConeTypeWeights(theNearbyTargetRGCindex);
        ratio = theCenterConeTypeWeights(2)/theCenterConeTypeWeights(1);
        if (theNearbyMajorityConeType ~= theMajorityConeType) && (ratio < 0.7)
            notFound = false;
            theNearbyTargetRGCindexOfDifferentMajorityConeType = theNearbyTargetRGCindex;
        end
    end
    if (notFound)
        theNearbyTargetRGCindexOfDifferentMajorityConeType = theTargetRGCindex;
        theNearbyMajorityConeType = theMajorityConeType;
    end
end