function generateSamplingGrids(obj, visualizeSpatialSamplingGrids)

    centerConnectableConeTypes = obj.rfModelParams.centerConnectableConeTypes;

    % Find out the range of cones in the RF center
    allConesNumPooledByTheRFcenters = full(sum(obj.theRGCMosaic.rgcRFcenterConeConnectivityMatrix,1));
    conesNumPooledByTheRFcenters = unique(allConesNumPooledByTheRFcenters);
    fprintf('Cones/RF center for this mosaic: %d\n', conesNumPooledByTheRFcenters);

    % Generate grid of spatial positions x conesNumInRFcenter
    [multiFocalOpticalPositionGrid, multiFocalConesNumInRFcenterGrid] = ...
        meshgrid(1:size(obj.nominalSpatialSamplingGrid,1), conesNumPooledByTheRFcenters);
    
    % Allocate memory
    multifocalRTVFobjectsNum = numel(multiFocalConesNumInRFcenterGrid);
    obj.opticalPositionGrid = zeros(multifocalRTVFobjectsNum, 2);
    obj.conesNumPooledByTheRFcenterGrid = zeros(multifocalRTVFobjectsNum,1);
    obj.visualSTFSurroundToCenterRcRatioGrid = zeros(multifocalRTVFobjectsNum, 1);
    obj.visualSTFSurroundToCenterIntegratedSensitivityRatioGrid = zeros(multifocalRTVFobjectsNum,1);


    for iMultifocalRTVFobjIndex = 1:multifocalRTVFobjectsNum
        % Sampling position (within the  mosaic)
        eccPositionIndex = multiFocalOpticalPositionGrid(iMultifocalRTVFobjIndex);

        % Start with the nominal position
        opticalPositionDegs = obj.nominalSpatialSamplingGrid(eccPositionIndex,:);

        % Cones num in RF center
        conesNumPooled = multiFocalConesNumInRFcenterGrid(iMultifocalRTVFobjIndex);

        % Get indices of center cones for the RGC that is closest to the opticalPositionDegs
        indicesOfRGCsWithThisManyCenterCones = find(allConesNumPooledByTheRFcenters == conesNumPooled);
        [~,idx] = min(sum((bsxfun(@minus, obj.theRGCMosaic.rgcRFpositionsDegs(indicesOfRGCsWithThisManyCenterCones,:), opticalPositionDegs)).^2,2));
        theTargetRGCindex = indicesOfRGCsWithThisManyCenterCones(idx);

        % Update the opticalPositionDegs for this RTVFobj to reflect the actual position of theTargetRGCindex
        opticalPositionDegs = obj.theRGCMosaic.rgcRFpositionsDegs(theTargetRGCindex,:);

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
        obj.opticalPositionGrid(iMultifocalRTVFobjIndex,:) = opticalPositionDegs;
        obj.conesNumPooledByTheRFcenterGrid(iMultifocalRTVFobjIndex) = conesNumPooled;
        obj.visualSTFSurroundToCenterRcRatioGrid(iMultifocalRTVFobjIndex) = surroundToCenterRcRatio;
        obj.visualSTFSurroundToCenterIntegratedSensitivityRatioGrid(iMultifocalRTVFobjIndex) = scIntSensitivity;

        % Save the target RGC index for this RTVF
        obj.targetRGCindex(iMultifocalRTVFobjIndex) = theTargetRGCindex;

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
