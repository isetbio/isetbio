function MosaicConnectorFull

    recomputePhase1 = ~true;
    tmpDir = '/Volumes/SamsungT3/MATLAB/toolboxes/isetbio/calculators/mosaicConnector';
    
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
        connectivityRadiusDeg = 15;

        % Load data for the analyzed region
        [RGCRFPositionsMicrons, RGCRFSpacingsMicrons, conePositionsMicrons, desiredConesToRGCratios] = ...
            loadData(whichEye, mosaicFOVDegs, eccentricitySamplesNumCones, ...
            eccentricitySamplesNumRGC, maxMovementPercentileCones, ...
            maxMovementPercentileRGC, bestIterationForConeMosaic,  ...
            bestIterationForRGCMosaic, connectivityRadiusDeg);

        % Compute connection matrix between the 2 mosaics
        save(fullfile(tmpDir,'tmp2.mat'), 'RGCRFPositionsMicrons', 'conePositionsMicrons', 'RGCRFSpacingsMicrons', 'desiredConesToRGCratios');
    else
       
        
    % Options
    orphanRGCpolicy = 'remove' ; % valid options: {'remove', 'share input'}
        
    % Treshold (x mean spacing) for removing cones/rgcs that are too close
    thresholdFractionForMosaicIncosistencyCorrection = 0.5;
        
    horizEccDegs = [0 0];
    fovDegs = [0.5 0.2];
        
    % Specify center in microns
    roi.center = round(1000*WatsonRGCModel.rhoDegsToMMs(horizEccDegs));
    roi.size = round(1000*WatsonRGCModel.rhoDegsToMMs(fovDegs));
    roi.margin = 5;

        
        
    visualizeConnectionProcess = ~true;
    visualizeMosaic = ~true;
        
    % Instantiate a plotlab object
    plotlabOBJ = plotlab();

    % Apply the default plotlab recipe overriding 
    % the color order and the figure size
    figHeightInches = 12;
    plotlabOBJ.applyRecipe(...
        'renderer', 'painters', ... %'opengl', ...
        'axesBox', 'on', ...
        'colorOrder', [0 0 0; 1 0 0.5], ...
        'axesTickLength', [0.015 0.01]/2,...
        'axesFontSize', 22, ...
        'figureWidthInches', figHeightInches*(roi.size(1)-2*roi.margin)/(roi.size(2)-2*roi.margin), ...
        'figureHeightInches', figHeightInches);

    computeStep1 = ~true;
    if (computeStep1)
        load(fullfile(tmpDir,'tmp2.mat'), 'RGCRFPositionsMicrons', 'conePositionsMicrons', 'RGCRFSpacingsMicrons', 'desiredConesToRGCratios');
        % Step 1. Remove inconstencies within and across mosaics
        [conePositionsMicrons, coneSpacingsMicrons,...
         RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
         desiredConesToRGCratios] = improveMosaicStats(conePositionsMicrons, ...
                       RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
                       desiredConesToRGCratios, ...
                       thresholdFractionForMosaicIncosistencyCorrection, roi);

        save(fullfile(tmpDir,'tmp2_step1.mat'), 'RGCRFPositionsMicrons', 'conePositionsMicrons', ...
            'RGCRFSpacingsMicrons', 'coneSpacingsMicrons', 'desiredConesToRGCratios');
        return;
    else
        load(fullfile(tmpDir,'tmp2_step1.mat'), 'RGCRFPositionsMicrons', 'conePositionsMicrons', ...
            'RGCRFSpacingsMicrons', 'coneSpacingsMicrons', 'desiredConesToRGCratios');
    end
    
    
    %visualizeRGCmosaic(90,RGCRFPositionsMicrons, RGCRFSpacingsMicrons, roi, 'original', plotlabOBJ);
                   
    computeStep2 = ~true;
    if (computeStep2)    
        % Step 2. Assign types (L,M,S) in the cone mosaic
        tritanopicAreaDiameterMicrons = 0.25 * 300;
        relativeSconeSpacing = 2.7;  % This results to around 8-9% S-cones
        LtoMratio = 2.0;  % Valid range: [0 - Inf]
        coneTypes = assignConeTypes(conePositionsMicrons, coneSpacingsMicrons, ...
            tritanopicAreaDiameterMicrons, relativeSconeSpacing, LtoMratio, roi, visualizeMosaic, plotlabOBJ);
        save(fullfile(tmpDir,'tmp2_step2.mat'), 'coneTypes');
        return;
    else
        load(fullfile(tmpDir,'tmp2_step2.mat'), 'coneTypes');
    end
    
    
    computeStep3 = ~true;
    if (computeStep3)  
    % Step 3. Connect cone to the midget RGC mosaic
        [midgetRGCconnectionMatrix, ...
         RGCRFPositionsMicrons, ...
         RGCRFSpacingsMicrons] = computeConnectionMatrix(...
                RGCRFPositionsMicrons, conePositionsMicrons, ...
                RGCRFSpacingsMicrons, coneSpacingsMicrons, ...
                coneTypes, desiredConesToRGCratios, ...
                orphanRGCpolicy, visualizeConnectionProcess);
        save(fullfile(tmpDir,'tmp2_step3.mat'), 'midgetRGCconnectionMatrix', ...
                    'RGCRFPositionsMicrons', 'RGCRFSpacingsMicrons');
        return;
    else
        load(fullfile(tmpDir,'tmp2_step3.mat'), 'midgetRGCconnectionMatrix', ...
                    'RGCRFPositionsMicrons', 'RGCRFSpacingsMicrons');
    end
    
    zLevels = [0.3 1 ];
    whichLevelsToContour = [1];
    
    depictRGCTesselation = true;
    if (depictRGCTesselation)
        displayIDs = ~true;
        fitEllipse = true;
        visualizeRFs(zLevels, whichLevelsToContour, midgetRGCconnectionMatrix, conePositionsMicrons, ...
            RGCRFPositionsMicrons, coneSpacingsMicrons, coneTypes, roi, fitEllipse, displayIDs, plotlabOBJ);
    end
    
    %visualizeRGCmosaic(91,RGCRFPositionsMicrons, RGCRFSpacingsMicrons, roi, 'final', plotlabOBJ);
         
    % Compute RF sizes via ellipse fitting
    [semiAxes,rfCenters] = computeRFsizes(zLevels, whichLevelsToContour, midgetRGCconnectionMatrix, conePositionsMicrons, RGCRFPositionsMicrons, coneSpacingsMicrons, coneTypes, roi);

    hFig = figure(222);
    eccNeuronsMicrons = sqrt(sum(rfCenters.^2,2));
    eccNeuronsDegs = WatsonRGCModel.rhoMMsToDegs(eccNeuronsMicrons/1000);
    rfCenterRadiusDegs = WatsonRGCModel.rhoMMsToDegs(mean(semiAxes,2)/1000);
    plot(eccNeuronsDegs, rfCenterRadiusDegs, '.'); hold on;
    set(gca, 'XLim', [0 10], 'YLim', [0 0.3]);
    drawnow;

    fName = sprintf('SemiAxes');
    plotlabOBJ.exportFig(hFig, 'png', fName, fullfile(pwd(), 'exports'));

        
end

