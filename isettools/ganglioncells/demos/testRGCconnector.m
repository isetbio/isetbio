function testRGCconnector

    rng(1);
    theInputConeMosaic = cMosaic(...
        'sizeDegs', 2*[1.0 0.2], ...
        'eccentricityDegs', [2 0], ...
        'coneDensities', [0.6 0.3 0.1]);

    instantiationMode = 'default';
    instantiationMode = 'custom cone-to-RGC density';
    %instantiationMode = 'custom RGC position matrix';

    switch (instantiationMode)
        case 'default'
            % Default instantiation, using mRGC mosaic
            rc = RGCconnector(theInputConeMosaic);

        case 'custom cone-to-RGC density'
            % Set desired cone-RGC density
            coneToRGCDensityRatio = 5;

            % Instantiation with custom density regular hex RGC mosaic
            rc = RGCconnector(theInputConeMosaic, ...
                'coneToRGCDensityRatio', coneToRGCDensityRatio);

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
            rc = RGCconnector(theInputConeMosaic, 'RGCRFpositionsMicrons', testRGCpositionsMicrons);
    end % switch

end


function testRGCconnectorOLD
    % Generate the input cone mosaic
    theInputConeMosaic = cMosaic(...
        'sizeDegs', 2*[0.5 0.5], ...
        'eccentricityDegs', [2 0]+10, ...
        'coneDensities', [0.6 0.3 0.1]);

%     theInputConeMosaic = cMosaic(...
%         'sizeDegs', [0.5 0.5], ...
%         'eccentricityDegs', [1 1], ...
%         'coneDensities', [0.6 0.3 0.1]);

    % Initialize with a precomputed RGC lattice
    RGCRFposMicrons = RGCRFconnector.initializeWithPrecomputedLattice(...
        theInputConeMosaic);


    % Initialize with a perfect hex lattice
    %coneToRGCDensityRatio = size(theInputConeMosaic.coneRFpositionsMicrons,1)/size(RGCRFposMicronsPrecomputed,1);
    %RGCRFposMicrons = RGCRFconnector.initializeWithPerfectHexLattice(...
    %    theInputConeMosaic, coneToRGCDensityRatio);
    %RGCRFposMicronsPerfect = RGCRFposMicrons;

    %visualizeInitialPositions(theInputConeMosaic, RGCRFposMicronsPerfect, RGCRFposMicronsPrecomputed);

    % Compute localConeToRGCDensityRatio struct, which we pass around 
    % so we can compute the coneToRGCdensityRatio at
    sampleSizeMicrons = 3;
    localConeToRGCDensityRatioStruct = RGCRFconnector.localConeToRGCDensityRatioFunctionHandle(...
            RGCRFposMicrons, theInputConeMosaic, sampleSizeMicrons);

    % Visualize the cone/rgc density map
    RGCRFconnector.visualizeLocalConeToRGCDensityMaps(...
        localConeToRGCDensityRatioStruct, ...
        RGCRFposMicrons, theInputConeMosaic);


    % Initial wiring, where each cone gets assigned to at least one RGC
    % If 2 or more cones get assigned to the same RGC, there can be a bias towards
    % minimizing spatial variance or minimizing chromatic variance.
    % This is controlled by the value of chromaticSpatialVarianceTradeoff, w,  
    % which is used to compute the cost for maintaning one RGC's cone inputs as follows:
    %
    % w-VALUE   COST FUNCTION                                           RESULT 
    % 1         proportional to spatial variance                        Minimize spatial variance, ignoring homogeneity of input cone types
    % 0         proportional to chromatic variance, e.g.:
    %             lConesNum/mConesNum, if lConesNum < mConesNum, or     Minimize chromatic variance, ignoring spatial position of input cones
    %             mConesNum/lConesNum, if mConesNum < lConesNum
    % 0<w<1     w * spatialVariance + (1-w)*chromatic variance          Minimize the combined chromatic+spatial variance

    wiringScenarios = {...
        'minimize spatial variance' ...
        'minimize spatial variance, sequential cone wiring' ...
        'minimize chromatic variance' ...
        'minimize chromatic variance, sequential cone wiring'};
    wiringScenarios = {wiringScenarios{3}};

    % Visualization parameters
    visualizationParams = struct(...
        'generateProgressionVideo', true, ...
        'exportStagesAsPDF', true, ...
        'generateSummaryFigure', true);

    for iScenario = 1:numel(wiringScenarios)
        doIt(wiringScenarios{iScenario}, theInputConeMosaic, RGCRFposMicrons, ...
            localConeToRGCDensityRatioStruct, visualizationParams);
    end
end




function doIt(wiringScenario, theInputConeMosaic, RGCRFposMicrons, ...
    localConeToRGCDensityRatioStruct, visualizationParams)

    switch wiringScenario
        case 'minimize spatial variance'
            chromaticSpatialVarianceTradeoff = 1;
            sequentialConeTypeWiring = false;

        case 'minimize spatial variance, sequential cone wiring'
            chromaticSpatialVarianceTradeoff = 1;
            sequentialConeTypeWiring = true;

        case 'minimize chromatic variance'
            chromaticSpatialVarianceTradeoff = 0;
            sequentialConeTypeWiring = false;

        case 'minimize chromatic variance, sequential cone wiring'
            chromaticSpatialVarianceTradeoff = 0;
            sequentialConeTypeWiring = true;
    end

    % Wiring parameters
    wiringParams = struct(...
        'chromaticSpatialVarianceTradeoff', chromaticSpatialVarianceTradeoff, ...
        'sequentialConeTypeWiring', sequentialConeTypeWiring, ...
        'maxNearbyRGCsNum', 7);

    % Wire mosaic
    [RGCRFinputs, RGCRFweights] = RGCRFconnector.wireInputConeMosaicToRGCRFs(...
        theInputConeMosaic, RGCRFposMicrons, ...
        localConeToRGCDensityRatioStruct, ...
        wiringParams,visualizationParams);
   

    if (visualizationParams.generateSummaryFigure)
    % Summary figure showing wired mosaic + stats
        hFig = figure(1000); clf;
        set(hFig, 'Color', [1 1 1], 'Position', [10 10 900 1100]);
        if (wiringParams.sequentialConeTypeWiring)
            set(hFig, 'Color', [1 1 1], 'Name', sprintf('Mosaic_W_%2.2f_SequentialConeWiring', wiringParams.chromaticSpatialVarianceTradeoff));
        else
            set(hFig, 'Color', [1 1 1], 'Name', sprintf('Mosaic_W_%2.2f_NonSequentialConeWiring', wiringParams.chromaticSpatialVarianceTradeoff));
        end
        axMosaic = subplot('Position', [0.15 0.35 0.8 0.62]);
        axStats{1} = subplot('Position', [0.07 0.05 0.42 0.2]);
        axStats{2} = subplot('Position', [0.57 0.05 0.42 0.2]);

        RGCRFconnector.visualizeWiringOfInputConeMosaicToRFcenters(...
            RGCRFinputs, RGCRFweights, ...
            theInputConeMosaic, ...
            'titleString', sprintf('wiring scenario: %s', wiringScenario), ...
            'superimposeConeInputWiring', true, ...
            'pauseAfterEachRGCisRendered', ~true, ...
            'figureHandle', hFig, ...
            'axesHandle', axMosaic);

        RGCRFconnector.reportConeInputStatistics(...
            wiringParams.chromaticSpatialVarianceTradeoff, ...
            RGCRFinputs, RGCRFweights, theInputConeMosaic, ...
            localConeToRGCDensityRatioStruct, ...
            'figureHandle', hFig, ...
            'axesHandles', axStats);

        NicePlot.exportFigToPDF(sprintf('Summary_%s.pdf', wiringScenario), hFig, 300);
    end

    
    
end

function visualizeInitialPositions()
    hFig = figure(1); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1400 600]);
    ax = subplot(1,2,1);
    [~, ~, XLims, YLims] = RGCRFconnector.plotRGCRFpos(RGCRFposMicronsPerfect,...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'inputConeMosaic', theInputConeMosaic, ...
        'thetaSamples', 30, ...
        'titleString', sprintf('perfect hex RGC lattice with corresponding density\n(RGF RFs num: %d)', size(RGCRFposMicronsPerfect,1)));
    ax = subplot(1,2,2);
    RGCRFconnector.plotRGCRFpos(RGCRFposMicronsPrecomputed,...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'inputConeMosaic', theInputConeMosaic, ...
        'XLims', XLims, 'YLims', YLims, ...
        'titleString', sprintf('mRGC RF lattice at corresponding ecc\n(RGF RFs num: %d)', size(RGCRFposMicronsPrecomputed,1)));
    
end
