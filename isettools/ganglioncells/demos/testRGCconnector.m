function testRGCconnector

    rng(1);

    eccentricity = 'very high';
    eccentricity = 'high';
    eccentricity = 'medium high';
    %eccentricity = 'medium';
    %eccentricity = 'medium low';
    %eccentricity = 'low';
    %eccentricity = 'very low';
    %eccentricity = 'foveal';

    switch (eccentricity)
        case 'very high'
            theInputConeMosaic = cMosaic(...
                'sizeDegs', 0.5*15*[0.6 0.4], ...
                'eccentricityDegs', [2 1]+20, ...
                'coneDensities', [0.6 0.3 0.1]);
        case 'high'
            theInputConeMosaic = cMosaic(...
                'sizeDegs', 0.5*10*[0.6 0.4], ...
                'eccentricityDegs', [2 1]+10, ...
                'coneDensities', [0.6 0.3 0.1]);

        case 'medium high'
            theInputConeMosaic = cMosaic(...
                'sizeDegs', 0.75*5*[0.6 0.4], ...
                'eccentricityDegs', [3 2]+5, ...
                'coneDensities', [0.6 0.3 0.1]);
        case 'medium'
            theInputConeMosaic = cMosaic(...
                'sizeDegs', 0.75*3*[0.6 0.4], ...
                'eccentricityDegs', [3 2]+2, ...
                'coneDensities', [0.6 0.3 0.1]);

        case 'medium low'
            theInputConeMosaic = cMosaic(...
                'sizeDegs', 0.75*1.5*[0.6 0.4], ...
                'eccentricityDegs', [1.5 1], ...
                'coneDensities', [0.6 0.3 0.1]);

        case 'low'
            theInputConeMosaic = cMosaic(...
                'sizeDegs', 0.75*1.5*[0.6 0.4], ...
                'eccentricityDegs', [1 0.6], ...
                'coneDensities', [0.6 0.3 0.1]);

        case 'very low'
            theInputConeMosaic = cMosaic(...
                'sizeDegs', 0.75*[0.6 0.4], ...
                'eccentricityDegs', [0.3 0.2], ...
                'coneDensities', [0.6 0.3 0.1]);

        case 'foveal'
            theInputConeMosaic = cMosaic(...
                'sizeDegs', 0.75*0.5*[0.6 0.4], ...
                'eccentricityDegs', [0 0], ...
                'coneDensities', [0.6 0.3 0.1]);

        otherwise
            error('Unknown eccentricity');
    end





    instantiationMode = 'default';
    %instantiationMode = 'custom cone-to-RGC density';
    %instantiationMode = 'custom RGC position matrix';

    % [0: minimize chromatic variance, 1: minimize spatial variance]
    chromaticSpatialVarianceTradeoff = 0.;
    maxNeighborNormDistance = 1.0;
    maxSwapPassesNum = 10;
    maxNumberOfConesToSwap = 6;
    visualizeIntermediateConnectivityStages = true;
    

    tic

    switch (instantiationMode)
        case 'default'
            % Default instantiation, using mRGC mosaic
            rc = RGCconnector(theInputConeMosaic, ...
                    'chromaticSpatialVarianceTradeoff', chromaticSpatialVarianceTradeoff, ...
                    'maxNeighborNormDistance', maxNeighborNormDistance, ...
                    'maxSwapPassesNum', maxSwapPassesNum, ...
                    'visualizeIntermediateConnectivityStages', visualizeIntermediateConnectivityStages);

        case 'custom cone-to-RGC density'
            % Set desired cone-RGC density
            coneToRGCDensityRatio = 5;

            % Instantiation with custom density regular hex RGC mosaic
            rc = RGCconnector(theInputConeMosaic, ...
                'coneToRGCDensityRatio', coneToRGCDensityRatio, ...
                'chromaticSpatialVarianceTradeoff', chromaticSpatialVarianceTradeoff, ...
                'maxNeighborNormDistance', maxNeighborNormDistance, ...
                'maxSwapPassesNum', maxSwapPassesNum, ...
                'visualizeIntermediateConnectivityStages', visualizeIntermediateConnectivityStages);

        case 'custom RGC position matrix'
            % Generate test RGC positions
            center = theInputConeMosaic.eccentricityMicrons;
            testRGCpositionsMicrons = [];
            for iCell = 0:6
                if (iCell == 0)
                    radius = 0;
                else
                    radius = 20;
                end
                testRGCpositionsMicrons(size(testRGCpositionsMicrons,1)+1,:) = center + radius*[cosd(iCell*60) sind(iCell*60)];
            end
            for iCell = 1:6
                radius = 50;
                testRGCpositionsMicrons(size(testRGCpositionsMicrons,1)+1,:) = center + radius*[cosd(iCell*60) sind(iCell*60)];
            end
    
            % Instantiation with custom RGC lattice positions
            rc = RGCconnector(theInputConeMosaic, ...
                'RGCRFpositionsMicrons', testRGCpositionsMicrons, ...
                'chromaticSpatialVarianceTradeoff', chromaticSpatialVarianceTradeoff, ...
                'maxNeighborNormDistance', maxNeighborNormDistance, ...
                'maxSwapPassesNum', maxSwapPassesNum, ...
                'visualizeIntermediateConnectivityStages', visualizeIntermediateConnectivityStages);
    end % switch

    toc

    % Export final connectivity
    hFig = rc.visualizeCurrentConnectivityState(9999);

    pfdFileName = sprintf('Ecc_%s_ChromaSpatialVarianceTradeoff_%2.2f.pdf',eccentricity, chromaticSpatialVarianceTradeoff);
    NicePlot.exportFigToPDF(pfdFileName, hFig, 300);
end

