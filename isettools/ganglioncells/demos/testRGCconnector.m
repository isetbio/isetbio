function testRGCconnector

    rng(1);

    %eccentricity = 'very high';
    eccentricity = 'high';
    %eccentricity = 'medium high';
    %eccentricity = 'medium';   % DONE
    %eccentricity = 'medium low';
    %eccentricity = 'low';
    %eccentricity = 'very low';
    %eccentricity = 'foveal';

    switch (eccentricity)
        case 'very high'
            theInputConeMosaic = cMosaic(...
                'sizeDegs', 0.425*15*[0.6 0.4], ...
                'eccentricityDegs', [17 18], ...
                'coneDensities', [0.6 0.3 0.1]);
        case 'high'
            theInputConeMosaic = cMosaic(...
                'sizeDegs', 0.5*10*[0.6 0.4], ...
                'eccentricityDegs', [12 11], ...
                'coneDensities', [0.6 0.3 0.1]);

        case 'medium high'
            theInputConeMosaic = cMosaic(...
                'sizeDegs', 0.7*5*[0.6 0.4], ...
                'eccentricityDegs', [8 7], ...
                'coneDensities', [0.6 0.3 0.1]);
        case 'medium'
            theInputConeMosaic = cMosaic(...
                'sizeDegs', 0.75*3*[0.6 0.4], ...
                'eccentricityDegs', [5 4], ...
                'coneDensities', [0.6 0.3 0.1]);

        case 'medium low'
            theInputConeMosaic = cMosaic(...
                'sizeDegs', 0.75*1.5*[0.6 0.4], ...
                'eccentricityDegs', [1.5 1], ...
                'coneDensities', [0.6 0.3 0.1]);

        case 'low'
            theInputConeMosaic = cMosaic(...
                'sizeDegs', 0.7*1.5*[0.6 0.4], ...
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

    rfOverlapFactor = 0.0;

    % Cone swapping phase params
    maxPassesNum = 30;
    maxNumberOfConesToSwap = 10;

    visualizeIntermediateConnectivityStages = ~true;
    
    tic

    % [0: minimize chromatic variance, 1: minimize spatial variance]
    wList = [0.0 0.2 0.3 0.4 0.5 0.7 1.0];
    searchRadiiList = [1.0];
    
    loadPreviouslyGeneratedRGCconnector = ~true;

    for iTradeOffIndex = 1:numel(wList)
        for iSearchIndex = 1:numel(searchRadiiList)

            close all
            chromaticSpatialVarianceTradeoff = wList(iTradeOffIndex);
            maxNeighborNormDistance = searchRadiiList(iSearchIndex);

            if (loadPreviouslyGeneratedRGCconnector)

                pfdFileName = sprintf('Ecc_%s_MaxNeighborDist_%2.2f_ChromaSpatialVarianceTradeoff_%2.2f.pdf',eccentricity, maxNeighborNormDistance, chromaticSpatialVarianceTradeoff);
                theRGCconnectorFileName = strrep(pfdFileName, '.pdf', '.mat');
                load(theRGCconnectorFileName, 'theRGCconnectorOBJ');
                fprintf('\nLoaded previously generated @RGCconnectorOBJ from %s\n', theRGCconnectorFileName);
        
                % Apply overlap factor
                rfOverlapFactor = 0.35;
                theRGCconnectorOBJ.expandRFsToOverlappingCones(...
                    'rfOverlapFactor', rfOverlapFactor ...
                    );

                % Visualize RFs with non-overlapping cones
                rgcsNum = size(theRGCconnectorOBJ.coneConnectivityMatrix,2);
                for iRGC = 1:rgcsNum
                    nonOverlappingConeIndices = find(squeeze(theRGCconnectorOBJ.coneConnectivityMatrix(:, iRGC))>0);
                    if (numel(nonOverlappingConeIndices)>1)
                        theRGCconnectorOBJ.visualizeConePoolingWithinRFcenter(iRGC);
                        pause
                    end
                end

                pfdFileName = strrep(pfdFileName, '.pdf', '');
                pfdFileName = sprintf('%s_rfOverlap_%2.2f.pdf', pfdFileName, rfOverlapFactor);
                continue;
            end

            switch (instantiationMode)
                case 'default'
                    % Default instantiation, using mRGC mosaic
                    theRGCconnectorOBJ = RGCconnector(theInputConeMosaic, ...
                            'rfOverlapFactor', rfOverlapFactor, ...
                            'chromaticSpatialVarianceTradeoff', chromaticSpatialVarianceTradeoff, ...
                            'maxNeighborNormDistance', maxNeighborNormDistance, ...
                            'maxPassesNum', maxPassesNum, ...
                            'maxNumberOfConesToSwap', maxNumberOfConesToSwap, ...
                            'visualizeIntermediateConnectivityStages', visualizeIntermediateConnectivityStages);
        
                case 'custom cone-to-RGC density'
                    % Set desired cone-RGC density
                    coneToRGCDensityRatio = 5;
        
                    % Instantiation with custom density regular hex RGC mosaic
                    theRGCconnectorOBJ = RGCconnector(theInputConeMosaic, ...
                        'coneToRGCDensityRatio', coneToRGCDensityRatio, ...
                        'rfOverlapFactor', rfOverlapFactor, ...
                        'chromaticSpatialVarianceTradeoff', chromaticSpatialVarianceTradeoff, ...
                        'maxNeighborNormDistance', maxNeighborNormDistance, ...
                        'maxPassesNum', maxPassesNum, ...
                        'maxNumberOfConesToSwap', maxNumberOfConesToSwap, ...
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
                    theRGCconnectorOBJ = RGCconnector(theInputConeMosaic, ...
                        'RGCRFpositionsMicrons', testRGCpositionsMicrons, ...
                        'rfOverlapFactor', rfOverlapFactor, ...
                        'chromaticSpatialVarianceTradeoff', chromaticSpatialVarianceTradeoff, ...
                        'maxNeighborNormDistance', maxNeighborNormDistance, ...
                        'maxPassesNum', maxPassesNum, ...
                        'maxNumberOfConesToSwap', maxNumberOfConesToSwap, ...
                        'visualizeIntermediateConnectivityStages', visualizeIntermediateConnectivityStages);
            end % switch
        
            toc
        
            % Export final connectivity
            hFig = theRGCconnectorOBJ.visualizeCurrentConnectivityState(9999);
            drawnow
           
            pfdFileName = sprintf('Ecc_%s_MaxNeighborDist_%2.2f_ChromaSpatialVarianceTradeoff_%2.2f.pdf',eccentricity, maxNeighborNormDistance, chromaticSpatialVarianceTradeoff);
            NicePlot.exportFigToPDF(pfdFileName, hFig, 300);
            
            % Export the generated @RGCconnector object
            theRGCconnectorFileName = strrep(pfdFileName, '.pdf', '.mat');
            save(theRGCconnectorFileName, 'theRGCconnectorOBJ', '-v7.3');
            fprintf('\nGenerated @RGCconnectorOBJ saved in %s\n', theRGCconnectorFileName);
        
        
        end % iSearchIndex
    end % iTradeOff
end

