function generateSamplingGrids(obj, minSpatialSamplingDegs, visualizeSpatialSamplingGrids)

    centerConnectableConeTypes = obj.rfModelParams.centerConnectableConeTypes;

    % Find out the range of cones in the RF center
    allConesNumPooledByTheRFcenters = full(sum(obj.theRGCMosaic.rgcRFcenterConeConnectivityMatrix,1));
    conesNumPooledByTheRFcenters = unique(allConesNumPooledByTheRFcenters);
    fprintf('Cones/RF center for this mosaic: %d\n', conesNumPooledByTheRFcenters);

    % Generate grid of spatial positions x conesNumInRFcenter
    [multiFocalOpticalPositionGrid, multiFocalConesNumInRFcenterGrid] = ...
        meshgrid(1:size(obj.nominalSpatialSamplingGrid,1), conesNumPooledByTheRFcenters);
    
    % Allocate memory
    maxmultifocalRTVFobjectsNum = numel(multiFocalConesNumInRFcenterGrid);
    spatialPositionsVisitedForConesNumPooled = cell(1, max(conesNumPooledByTheRFcenters));

    obj.opticalPositionGrid = [];
    obj.conesNumPooledByTheRFcenterGrid = [];
    obj.visualSTFSurroundToCenterRcRatioGrid = [];
    obj.visualSTFSurroundToCenterIntegratedSensitivityRatioGrid = [];
    obj.targetRGCindex = [];

    for iMultifocalRTVFobjIndex = 1:maxmultifocalRTVFobjectsNum

        % Cones num in RF center
        conesNumPooled = multiFocalConesNumInRFcenterGrid(iMultifocalRTVFobjIndex);

        % Sampling position (within the  mosaic)
        eccPositionIndex = multiFocalOpticalPositionGrid(iMultifocalRTVFobjIndex);

        % Start with the nominal position
        opticalPositionDegs = obj.nominalSpatialSamplingGrid(eccPositionIndex,:);

        % Get indices of center cones for the RGC that is closest to the opticalPositionDegs
        indicesOfRGCsWithThisManyCenterCones = find(allConesNumPooledByTheRFcenters == conesNumPooled);
        [~,idx] = min(sum((bsxfun(@minus, obj.theRGCMosaic.rgcRFpositionsDegs(indicesOfRGCsWithThisManyCenterCones,:), opticalPositionDegs)).^2,2));
        theTargetRGCindex = indicesOfRGCsWithThisManyCenterCones(idx);

        % Update the opticalPositionDegs for this RTVFobj to reflect the actual position of theTargetRGCindex
        opticalPositionDegs = obj.theRGCMosaic.rgcRFpositionsDegs(theTargetRGCindex,:);

        % Check to see if we visited a nearby position
        alreadyVisitedPositions = spatialPositionsVisitedForConesNumPooled{conesNumPooled};
        if (~isempty(alreadyVisitedPositions))
            mindistance = min(sqrt(sum((bsxfun(@minus, alreadyVisitedPositions, [opticalPositionDegs(1) opticalPositionDegs(2)])).^2,2)));
            if (mindistance < minSpatialSamplingDegs)
                continue;
            end
        end
       
        % We havent visited, so add it
        alreadyVisitedPositions = cat(1, alreadyVisitedPositions, [opticalPositionDegs(1) opticalPositionDegs(2)]);
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
        temporalEquivalentEccDegs = obj.theRGCMosaic.temporalEquivalentEccentricityForEccentricity(opticalPositionDegs);
        radialTemporalEquivalentEccDegs = sqrt(sum(temporalEquivalentEccDegs.^2,2));
        scIntSensitivity = RGCmodels.CronerKaplan.constants.surroundToCenterIntegratedSensitivityRatioFromEccDegsForPcells(radialTemporalEquivalentEccDegs);
       
        % The 4 grids
        obj.opticalPositionGrid = ...
            cat(1, obj.opticalPositionGrid, opticalPositionDegs);

        obj.conesNumPooledByTheRFcenterGrid = ...
            cat(1, obj.conesNumPooledByTheRFcenterGrid,conesNumPooled);

        obj.visualSTFSurroundToCenterRcRatioGrid  = ...
            cat(1, obj.visualSTFSurroundToCenterRcRatioGrid, surroundToCenterRcRatio);

        obj.visualSTFSurroundToCenterIntegratedSensitivityRatioGrid = ...
            cat(1, obj.visualSTFSurroundToCenterIntegratedSensitivityRatioGrid, scIntSensitivity);

         % Save the target RGC index for this RTVF
         obj.targetRGCindex = cat(1, obj.targetRGCindex, theTargetRGCindex);
    end % for iMultifocalRTVFobjIndex 

    if (visualizeSpatialSamplingGrids)
        hFig = figure(2); clf;
        set(hFig, 'Color', [1 1 1]);

        for iConesNumPooled = 1:numel(conesNumPooledByTheRFcenters)
            ax = subplot(1,numel(conesNumPooledByTheRFcenters), iConesNumPooled);

            % Determine spatial grid coords for this # of center cones
            conesNumPooled = conesNumPooledByTheRFcenters(iConesNumPooled);
            idx = find(obj.conesNumPooledByTheRFcenterGrid == conesNumPooled);
            thisSpatialGrid = obj.opticalPositionGrid(idx,:);

            obj.plotSpatialSamplingGrid(ax, thisSpatialGrid, sprintf('spatial sampling grid for %d center cones', conesNumPooled));
        end

    end
end
