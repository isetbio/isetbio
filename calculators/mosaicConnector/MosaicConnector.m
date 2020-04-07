function MosaicConnector

    recomputePhase1 = ~true;
    
    if (recomputePhase1)
        % Select mosaics to load
        whichEye = 'right';
        mosaicFOVDegs = 15;
        eccentricitySamplesNumCones = 32;  
        eccentricitySamplesNumRGC = 32; 
        maxMovementPercentileCones = 20;
        maxMovementPercentileRGC = 20;
        bestIterationForConeMosaic = Inf;
        bestIterationForRGCMosaic = 95;

        % Connect mosaics only within a central region to save compute time
        connectivityRadiusDeg = 12;

        % Load data for the analyzed region
        [RGCRFPositionsMicrons, RGCRFSpacingsMicrons, conePositionsMicrons, desiredConesToRGCratios] = ...
            loadData(whichEye, mosaicFOVDegs, eccentricitySamplesNumCones, ...
            eccentricitySamplesNumRGC, maxMovementPercentileCones, ...
            maxMovementPercentileRGC, bestIterationForConeMosaic,  ...
            bestIterationForRGCMosaic, connectivityRadiusDeg);

        % Compute connection matrix between the 2 mosaics
        save('tmp.mat', 'RGCRFPositionsMicrons', 'conePositionsMicrons', 'RGCRFSpacingsMicrons', 'desiredConesToRGCratios');
    else
        load('tmp.mat', 'RGCRFPositionsMicrons', 'conePositionsMicrons', 'RGCRFSpacingsMicrons', 'desiredConesToRGCratios');
        
        % Define region of interest to work on
        roi.center = [1900 0];
        roi.size = [100 50];
        
        % Instantiate a plotlab object
        plotlabOBJ = plotlab();

        % Apply the default plotlab recipe overriding 
        % the color order and the figure size
        figHeightInches = 10;
        plotlabOBJ.applyRecipe(...
            'colorOrder', [0 0 0; 1 0 0.5], ...
            'figureWidthInches', figHeightInches*roi.size(1)/roi.size(2), ...
            'figureHeightInches', figHeightInches);
    
        % Some times, elements are too close. Remove them
        thresholdSeparationMicronsForRemovingUnitsFromMosaic = 0.0;
    
        [connectivityMatrix, ...
         conePositionsMicrons, ...
         RGCRFPositionsMicrons, ...
         coneSpacingsMicrons] = computeConnectionMatrix(...
                RGCRFPositionsMicrons, conePositionsMicrons, ...
                RGCRFSpacingsMicrons, desiredConesToRGCratios, ...
                roi, thresholdSeparationMicronsForRemovingUnitsFromMosaic);

        visualizeRFs(connectivityMatrix, conePositionsMicrons, ...
            RGCRFPositionsMicrons, coneSpacingsMicrons, roi)
    end
    
end

