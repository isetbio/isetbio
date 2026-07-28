function dataOut = componentsForRFcenterConnectivity(...
    theInputConeMosaic, extraDegsToSupportRGCSurrounds, ...
    rfCenterConnectivityParams, theIntermediateConnectivityStagesMetaDataFileName)

    % Generate metaData struct for sourceLatticeStruct
    metaDataStruct.coneTypes = theInputConeMosaic.coneTypes;
    metaDataStruct.coneTypeIDs = [theInputConeMosaic.LCONE_ID    theInputConeMosaic.MCONE_ID    theInputConeMosaic.SCONE_ID];
    metaDataStruct.coneColors  = [theInputConeMosaic.lConeColor; theInputConeMosaic.mConeColor; theInputConeMosaic.sConeColor];
    metaDataStruct.midgetRGCSurroundRadiusMicronsAtMaxEccentricity = 1e3 * theInputConeMosaic.customDegsToMMsConversionFunction(extraDegsToSupportRGCSurrounds/2);
    
   % Generate sourceLatticeStruct for the @coneToMidgetRGCConnector
    sourceLatticeStruct = struct(...
            'name', 'cone RFs', ...
            'DegsToMMsConversionFunction', theInputConeMosaic.customDegsToMMsConversionFunction, ...
            'MMsToDegsConversionFunction', theInputConeMosaic.customMMsToDegsConversionFunction, ...
            'RFpositionsMicrons', theInputConeMosaic.coneRFpositionsMicrons, ...
            'metaData', metaDataStruct ...
    );
    
    % Import mRGC RF positions (destination) 
    mRGCRFposMicrons = retinalattice.import.finalMRGCPositions(...
         theInputConeMosaic.sourceLatticeSizeDegs, ...
         mean(sourceLatticeStruct.RFpositionsMicrons,1), ... 
         max(max(sourceLatticeStruct.RFpositionsMicrons,[], 1)-min(sourceLatticeStruct.RFpositionsMicrons,[], 1)), ...
         theInputConeMosaic.whichEye, ...
         theInputConeMosaic.customDegsToMMsConversionFunction);

    % Generate the destination lattice struct
    destinationLatticeStruct = struct(...
        'name', 'mRGC RFs', ...
        'DegsToMMsConversionFunction', theInputConeMosaic.customDegsToMMsConversionFunction, ...
        'MMsToDegsConversionFunction', theInputConeMosaic.customMMsToDegsConversionFunction, ...
        'RFpositionsMicrons', mRGCRFposMicrons ...
        );

    % Cone indices to be connected
    coneIndicesToBeConnected = [];
    if (ismember(cMosaic.LCONE_ID, rfCenterConnectivityParams.coneTypesToBeConnected))
        fprintf('L-cones allowed to be connected to the RF center\n');
        coneIndicesToBeConnected = cat(1, coneIndicesToBeConnected, theInputConeMosaic.lConeIndices(:));
    end
    if (ismember(cMosaic.MCONE_ID, rfCenterConnectivityParams.coneTypesToBeConnected))
        fprintf('M-cones allowed to be connected to the RF center\n');
        coneIndicesToBeConnected = cat(1, coneIndicesToBeConnected, theInputConeMosaic.mConeIndices(:));
    end
    if (ismember(cMosaic.SCONE_ID, rfCenterConnectivityParams.coneTypesToBeConnected))
        fprintf('S-cones allowed to be connected to the RF center\n');
        coneIndicesToBeConnected = cat(1, coneIndicesToBeConnected, theInputConeMosaic.sConeIndices(:));
    end


    % Connect the source and destination lattices
    theMosaicConnectorOBJ = coneToMidgetRGCConnector(...
        sourceLatticeStruct, destinationLatticeStruct, ...
        'optimizationCenter', rfCenterConnectivityParams.optimizationCenter, ...
        'maxNeighborNormDistance', rfCenterConnectivityParams.maxNeighborNormDistance, ...
        'maxNeighborsNum', rfCenterConnectivityParams.maxNeighborsNum, ...
        'maxConeInputsPerRGCToConsiderTransferToNearbyRGCs', rfCenterConnectivityParams.maxConeInputsPerRGCToConsiderTransferToNearbyRGCs, ...
        'maxConeInputsPerRGCToConsiderSwappingWithNearbyRGCs', rfCenterConnectivityParams.maxConeInputsPerRGCToConsiderSwappingWithNearbyRGCs, ...
        'maxPassesNum', rfCenterConnectivityParams.maxPassesNum, ...
        'localSpacingFromCurrentCentroids', rfCenterConnectivityParams.localSpacingFromCurrentCentroids, ...
        'spatialChromaticUniformityTradeoff', rfCenterConnectivityParams.spatialChromaticUniformityTradeoff, ...
        'coneIndicesToBeConnected', coneIndicesToBeConnected, ...
        'saveIntermediateConnectivityStagesMetaData', rfCenterConnectivityParams.saveIntermediateStagesOfCenterConnectivityOptimization, ...
        'visualizeConnectivityAtIntermediateStages', rfCenterConnectivityParams.visualizeIntermediateStagesOfCenterConnectivityOptimization, ...
        'generateProgressVideo', ~true);


    dataOut.rfCenterConnectivityParams = rfCenterConnectivityParams;
    dataOut.inputConeMosaic = theInputConeMosaic;

    % Save the intermediate metaData structs
    if (rfCenterConnectivityParams.saveIntermediateStagesOfCenterConnectivityOptimization)
        intermediateMetaDataStructs = theMosaicConnectorOBJ.intermediateMetaDataStructs;
        save(theIntermediateConnectivityStagesMetaDataFileName, 'intermediateMetaDataStructs', '-v7.3');
        % Save the intermediate figure handles
        %dataOut.mosaicConnectorIntermediateFigureHandles = theMosaicConnectorOBJ.intermediateFigureHandles;
    end

    % Set the rf center connectivity matrix
    dataOut.rgcRFcenterConeConnectivityMatrix = theMosaicConnectorOBJ.connectivityMatrix;
                
    % Set the rf positions, in microns
    dataOut.rgcRFpositionsMicrons = theMosaicConnectorOBJ.destinationRFcentroidsFromInputs;

    % Set the rf spacings, in microns
    dataOut.rgcRFspacingsMicrons = theMosaicConnectorOBJ.destinationRFspacingsFromCentroids;
   
    % Set the rf positions in degrees
    dataOut.rgcRFpositionsDegs = theInputConeMosaic.customMMsToDegsConversionFunction(dataOut.rgcRFpositionsMicrons/1000);
    
    % Set the rf spacings in degrees
    RFeccentricityMicrons = sqrt(sum((dataOut.rgcRFpositionsMicrons).^2,2));
    dataOut.rgcRFspacingsDegs = reshape(...
        theInputConeMosaic.sizeMicronsToSizeDegreesForCmosaic(dataOut.rgcRFspacingsMicrons(:), RFeccentricityMicrons(:)), size(dataOut.rgcRFspacingsMicrons));

    % Crop RGCs on the border
    dataOut = cropRGCsOnTheBorder(dataOut);
end


function dataOut = cropRGCsOnTheBorder(dataOut)
    % Remove RGC RFs on the margins 
    indicesOfRGCRFsToKeep = findRGCRFsWithinBorders(dataOut.rgcRFpositionsMicrons, dataOut.rgcRFspacingsMicrons);
    
    % Update rf center connectivity matrix & rgcRFpositionsMicrons
    dataOut.rgcRFcenterConeConnectivityMatrix = dataOut.rgcRFcenterConeConnectivityMatrix(:,indicesOfRGCRFsToKeep);

    dataOut.rgcRFpositionsMicrons = dataOut.rgcRFpositionsMicrons(indicesOfRGCRFsToKeep,:);
    dataOut.rgcRFspacingsMicrons = dataOut.rgcRFspacingsMicrons(indicesOfRGCRFsToKeep);
    dataOut.rgcRFpositionsDegs = dataOut.rgcRFpositionsDegs(indicesOfRGCRFsToKeep,:);
    dataOut.rgcRFspacingsDegs = dataOut.rgcRFspacingsDegs(indicesOfRGCRFsToKeep);

    dataOut.rgcsNum = numel(dataOut.rgcRFspacingsDegs);
end

function indicesOfRFsToKeep = findRGCRFsWithinBorders(rgcRFpositions, rgcRFspacings)

    idx = find(~isinf(rgcRFspacings));
    rgcRFpositions = rgcRFpositions(idx,:);
    rgcRFspacings = rgcRFspacings(idx);

    [minXY, idx] = min(rgcRFpositions,[],1);
    maxSpacing = max(rgcRFspacings(idx));
    [maxXY, idx] = max(rgcRFpositions,[],1);
    maxSpacing = max([maxSpacing max(rgcRFspacings(idx))]);

    xLims(1) = minXY(1)+maxSpacing;
    xLims(2) = maxXY(1)-maxSpacing;
    yLims(1) = minXY(2)+maxSpacing;
    yLims(2) = maxXY(2)-maxSpacing;
    indicesOfRFsToKeep= find(...
            (rgcRFpositions(:,1)>=xLims(1)) & ...
            (rgcRFpositions(:,1)<=xLims(2)) & ...
            (rgcRFpositions(:,2)>=yLims(1)) & ...
            (rgcRFpositions(:,2)<=yLims(2)) ...
            );
end