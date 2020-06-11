function computeWeightedConeInputsToRGCCenterSurroundSubregions(synthesizedRFParams, dataSet, ...
    conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
    midgetRGCconnectionMatrixCenter, ...
    plotlabOBJ, pdfFileName, exportsDir)

    % Get data
    switch (dataSet)
        case 'visual'
            rfParams = synthesizedRFParams.visual;
        case 'retinal'
            rfParams = synthesizedRFParams.retinal;
        otherwise
            error('Unknown dataSet: ''%s''.', dataSet)
    end
    
    % Ecc and indices of RGCs in the analyzed patch
    eccDegs = synthesizedRFParams.eccDegs;
    rgcIndices = synthesizedRFParams.rgcIndices;
    
    % Analyzed patch
    patchEccDegs = synthesizedRFParams.patchEccDegs;
    patchSizeDegs = synthesizedRFParams.patchSizeDegs;

    hFig = figure(3);
    
    % Allocate sparse matrix to save the connections of cones to RGC surrounds
    conesNum = size(midgetRGCconnectionMatrixCenter,1);
    rgcsNum = size(midgetRGCconnectionMatrixCenter,2);
    maxNumberOfConnections = rgcsNum*350;
    midgetRGCconnectionMatrixSurround = spalloc(conesNum, rgcsNum, maxNumberOfConnections);
    
     maxGainDivisor = 10;
     
    for iRGC = ceil(97*numel(rgcIndices)/100):numel(rgcIndices)

        % Get index of RGC in full mosaic
        rgcIndex = synthesizedRFParams.rgcIndices(iRGC);
        
        % Get position of RGC in microns
        rgcPositionMicrons = synthesizedRFParams.rgcCenterPositionMicrons(iRGC,:);
        
        % Get radius of RGC center in microns
        rgcCenterRadiusDegs = rfParams.centerRadiiDegs(iRGC);
        rgcCenterRadiusMicrons = convertSizeDegsToSizeMicrons(rgcCenterRadiusDegs, eccDegs(iRGC));
        
        % Get peak sensitivity of RGC center
        rgcCenterPeakSensivity = rfParams.centerPeakSensitivities(iRGC);
        
        % Get radius of RGC surround in microns
        rgcSurroundRadiusDegs = rfParams.surroundRadiiDegs(iRGC);
        rgcSurroundRadiusMicrons = convertSizeDegsToSizeMicrons(rgcSurroundRadiusDegs, eccDegs(iRGC));
        
        % Get peak sensitivity of RGC surround
        rgcSurroundPeakSensivity = rfParams.surroundPeakSensitivities(iRGC);
        
        % Retrieve cone indices connected to the RGC surround
        weights = gaussianConeWeights(conePositionsMicrons, rgcPositionMicrons, rgcSurroundRadiusMicrons);
        coneIndicesConnectedToSurround = find(weights >= exp(-3));
        midgetRGCconnectionMatrixSurround(coneIndicesConnectedToSurround, rgcIndex) = weights(coneIndicesConnectedToSurround) * rgcSurroundPeakSensivity;
        
        % Retrieve cone indices connected to the RGC center
        connectivityVector = full(squeeze(midgetRGCconnectionMatrixCenter(:, rgcIndex)));
        coneIndicesConnectedToCenter = find(connectivityVector>0);
        weights = gaussianConeWeights(conePositionsMicrons(coneIndicesConnectedToCenter,:), rgcPositionMicrons, rgcCenterRadiusMicrons);
        midgetRGCconnectionMatrixCenter(coneIndicesConnectedToCenter, rgcIndex) = weights * rgcCenterPeakSensivity;
        
        % center weights
        centerWeights = squeeze(full(midgetRGCconnectionMatrixCenter(coneIndicesConnectedToCenter, rgcIndex)));
        
        % surround weights
        surroundWeights = squeeze(full(midgetRGCconnectionMatrixSurround(coneIndicesConnectedToSurround, rgcIndex)));
        
        % Only visualize weights up to max(center)/300. After that saturate
        % to max
        maxWeightVisualized = max(centerWeights) / maxGainDivisor; % max([max(centerWeights) max(surroundWeights)]);
        
        % Plot center weights
        ax = subplot(1,2,1);
        visualizeWeigtedConeInputsToRGCSubregion(ax, coneTypes(coneIndicesConnectedToCenter), ...
            conePositionsMicrons(coneIndicesConnectedToCenter,:), ...
            coneSpacingsMicrons(coneIndicesConnectedToCenter), ...
            centerWeights,  maxWeightVisualized, ...
            rgcPositionMicrons, rgcCenterRadiusMicrons, rgcSurroundRadiusMicrons*2);
        
        % Plot surround weights
        ax = subplot(1,2,2);
        visualizeWeigtedConeInputsToRGCSubregion(ax, coneTypes(coneIndicesConnectedToSurround), ...
            conePositionsMicrons(coneIndicesConnectedToSurround,:), ...
            coneSpacingsMicrons(coneIndicesConnectedToSurround), ...
            surroundWeights,  maxWeightVisualized, ...
            rgcPositionMicrons, rgcSurroundRadiusMicrons, rgcSurroundRadiusMicrons*2);
        
        drawnow;
        pause
    end % iRGC

    plotlabOBJ.exportFig(hFig, 'png', sprintf('%s__%s',pdfFileName, dataSet),exportsDir);
end


function weights = gaussianConeWeights(conePositionsMicrons, rgcPositionMicrons, rgcSubregionRadiusMicrons)
    weights = exp(-sum((bsxfun(@minus, conePositionsMicrons, rgcPositionMicrons)).^2,2)/rgcSubregionRadiusMicrons^2);
end

function sizeMicrons = convertSizeDegsToSizeMicrons(sizeDegs, eccDegs)
    sizeMicrons = ...
        1000*(WatsonRGCModel.rhoDegsToMMs(eccDegs + sizeDegs/2) - WatsonRGCModel.rhoDegsToMMs(eccDegs - sizeDegs/2));
end
