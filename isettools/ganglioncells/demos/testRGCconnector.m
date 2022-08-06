function testRGCconnector

    rng(1);

    eccentricity = 'very high';
    eccentricity = 'high';
    eccentricity = 'medium high';
    eccentricity = 'medium';   % DONE
    %eccentricity = 'medium low';
    %eccentricity = 'low';
    %eccentricity = 'very low';
    %eccentricity = 'foveal';

     

    switch (eccentricity)
        case 'very high'
            theInputConeMosaic = cMosaic(...
                'sourceLatticeSizeDegs', 60, ...
                'overlappingConeFractionForElimination', 0.5, ...
                'sizeDegs', 0.475*15*[0.6 0.4], ...
                'eccentricityDegs', [17 18], ...
                'coneDensities', [0.6 0.3 0.1]);
        case 'high'
            theInputConeMosaic = cMosaic(...
                'sourceLatticeSizeDegs', 60, ...
                'overlappingConeFractionForElimination', 0.5, ...
                'sizeDegs', 0.5*10*[0.6 0.4], ...
                'eccentricityDegs', [12 11], ...
                'coneDensities', [0.6 0.3 0.1]);

        case 'medium high'
            theInputConeMosaic = cMosaic(...
                'sourceLatticeSizeDegs', 60, ...
                'overlappingConeFractionForElimination', 0.5, ...
                'sizeDegs', 0.7*5*[0.6 0.4], ...
                'eccentricityDegs', [8 7], ...
                'coneDensities', [0.6 0.3 0.1]);
        case 'medium'
            theInputConeMosaic = cMosaic(...
                'sourceLatticeSizeDegs', 60, ...
                'overlappingConeFractionForElimination', 0.5, ...
                'sizeDegs', 0.75*3*[0.6 0.4], ...
                'eccentricityDegs', [5 4], ...
                'coneDensities', [0.6 0.3 0.1]);

        case 'medium low'
            theInputConeMosaic = cMosaic(...
                'sourceLatticeSizeDegs', 60, ...
                'overlappingConeFractionForElimination', 0.5, ...
                'sizeDegs', 0.75*1.5*[0.6 0.4], ...
                'eccentricityDegs', [1.5 1], ...
                'coneDensities', [0.6 0.3 0.1]);

        case 'low'
            theInputConeMosaic = cMosaic(...
                'sourceLatticeSizeDegs', 60, ...
                'overlappingConeFractionForElimination', 0.5, ...
                'sizeDegs', 0.7*1.5*[0.6 0.4], ...
                'eccentricityDegs', [1 0.6], ...
                'coneDensities', [0.6 0.3 0.1]);

        case 'very low'
            theInputConeMosaic = cMosaic(...
                'sourceLatticeSizeDegs', 60, ...
                'overlappingConeFractionForElimination', 0.5, ...
                'sizeDegs', 0.75*[0.6 0.4], ...
                'eccentricityDegs', [0.3 0.2], ...
                'coneDensities', [0.6 0.3 0.1]);

        case 'foveal'
            theInputConeMosaic = cMosaic(...
                'sourceLatticeSizeDegs', 60, ...
                'overlappingConeFractionForElimination', 0.5, ...
                'sizeDegs', 0.75*0.5*[0.6 0.4], ...
                'eccentricityDegs', [0 0], ...
                'coneDensities', [0.6 0.3 0.1]);

        otherwise
            error('Unknown eccentricity');
    end




    loadPreviouslyGeneratedRGCconnector = ~true;
    
    instantiationMode = 'default';
    instantiationMode = 'custom cone-to-RGC density';
    %instantiationMode = 'custom RGC position matrix';

    % Cone swapping phase params
    maxPassesNum = 30;
    maxNumberOfConesToSwap = 8;
    maxMeanConeInputsPerRGCToConsiderSwapping = 10;

    visualizeIntermediateConnectivityStages = true;
    
    tic

    % [0: minimize chromatic variance, 1: minimize spatial variance]
    wList = [0.0]; %[0.0 0.2 0.35 0.5 0.65 0.8];
    maxNeighborNormDistanceList = [1.5];
    rfCentroidOverlapPenaltyFactor = 3;
    
    RcToRGCseparationRatio = 1.5;

    for iTradeOffIndex = 1:numel(wList)
        for iSearchIndex = 1:numel(maxNeighborNormDistanceList)

            close all
            chromaticSpatialVarianceTradeoff = wList(iTradeOffIndex);
            maxNeighborNormDistance = maxNeighborNormDistanceList(iSearchIndex);

            if (loadPreviouslyGeneratedRGCconnector)

                pfdFileName = sprintf('Ecc_%s_MaxNeighborDist_%2.2f_ChromaSpatialVarianceTradeoff_%2.2f_centroidOverlapPenalty_%2.2f.pdf',eccentricity, maxNeighborNormDistance, chromaticSpatialVarianceTradeoff, rfCentroidOverlapPenaltyFactor);
           
                theRGCconnectorFileName = strrep(pfdFileName, '.pdf', '.mat');
                load(theRGCconnectorFileName, 'theRGCconnectorOBJ');
                fprintf('\nLoaded previously generated @RGCconnectorOBJ from %s\n', theRGCconnectorFileName);
        
                % Apply overlap factor
                theRGCconnectorOBJ.divergeConeOutputsToMultipleNearbyRGCs(...
                    'RcToRGCseparationRatio', RcToRGCseparationRatio ...
                    );

                pfdFileName = strrep(pfdFileName, '.pdf', '');
                
                visualizedConesNum = [];
                switch (eccentricity)
                    case  'very high'
                        visualizedFieldOfViewMicrons = 50;
                    case 'high'
                        visualizedFieldOfViewMicrons = 150;
                        visualizedConesNum = 600;
                    case 'medium high'
                        visualizedFieldOfViewMicrons = 100;
                        visualizedConesNum = 400;
                    case 'medium'
                        visualizedFieldOfViewMicrons = 50;
                        visualizedConesNum = 100;
                    case 'medium low'
                        visualizedFieldOfViewMicrons = 30;
                        visualizedConesNum = 100;
                    case 'low'
                        visualizedFieldOfViewMicrons = 30;
                        visualizedConesNum = 100;
                    case 'very low'
                        visualizedFieldOfViewMicrons = 30;
                        visualizedConesNum = 100;
                    case'foveal'
                        visualizedFieldOfViewMicrons = 30;
                        visualizedConesNum = 100;
                end
                
    
                % Visualize the center-most RGC RF
                ecc = sqrt(sum((bsxfun(@minus, theRGCconnectorOBJ.RGCRFcentroidsFromInputs, mean(theRGCconnectorOBJ.RGCRFcentroidsFromInputs,1))).^2,2));
                [~,theCenterMostRGCindex] = min(ecc);
                
                [hFig,visualizedNeighbors] = theRGCconnectorOBJ.visualizeConePoolingWithinNeighboringRGCs(theCenterMostRGCindex, ...
                            'visualizedFieldOfViewMicrons', visualizedFieldOfViewMicrons, ...
                            'visualizedConesNum', visualizedConesNum);
                
                visualize1Doverlap = true;
                if (visualize1Doverlap)
    
                    for iNeighbor = 1:visualizedNeighbors
                        hFig = theRGCconnectorOBJ.visualizeConePoolingWithinRFcenter(theCenterMostRGCindex, ...
                                'visualizedFieldOfViewMicrons', visualizedFieldOfViewMicrons, ...
                                'visualizedConesNum', visualizedConesNum, ...
                                'visualizedNeighbor', iNeighbor);

                        pdfFileNameFinal = sprintf('%s_RcToRGCseparationRatio_%2.2f_Neighbor_%d.pdf',...
                            pfdFileName, RcToRGCseparationRatio, iNeighbor);

                        NicePlot.exportFigToPDF(pdfFileNameFinal, hFig, 300);
                    end
                end
                
                % Visualize all RGCRFs
                generateVideoShowingAllRGCs = false;
                if (generateVideoShowingAllRGCs)
                    movFileName = sprintf('%s_RFoverlapAnalysis', strrep(pdfFileName, '.pdf', ''));
                    videoOBJ = VideoWriter(movFileName, 'MPEG-4');
                    videoOBJ.FrameRate = 30;
                    videoOBJ.Quality = 100;
                    videoOBJ.open();
                    
                    rgcsNum = size(theRGCconnectorOBJ.coneConnectivityMatrix,2);
                    for iRGC = 1:rgcsNum
                        hFig = theRGCconnectorOBJ.visualizeConePoolingWithinRFcenter(iRGC, ...
                            'visualizedFieldOfViewMicrons', visualizedFieldOfViewMicrons, ...
                            'visualizedConesNum', visualizedConesNum);
                        drawnow;
                        videoOBJ.writeVideo(getframe(hFig));
                    end

                    videoOBJ.close();
                end
                
                continue;
            end

            % Generate new connected RGC mosaic
            switch (instantiationMode)
                case 'default'
                    % Default instantiation, using mRGC mosaic
                    theRGCconnectorOBJ = RGCconnector(theInputConeMosaic, ...
                            'chromaticSpatialVarianceTradeoff', chromaticSpatialVarianceTradeoff, ...
                            'rfCentroidOverlapPenaltyFactor', rfCentroidOverlapPenaltyFactor, ...
                            'maxNeighborNormDistance', maxNeighborNormDistance, ...
                            'maxPassesNum', maxPassesNum, ...
                            'maxMeanConeInputsPerRGCToConsiderSwapping', maxMeanConeInputsPerRGCToConsiderSwapping, ...
                            'maxNumberOfConesToSwap', maxNumberOfConesToSwap, ...
                            'visualizeIntermediateConnectivityStages', visualizeIntermediateConnectivityStages);
        
                case 'custom cone-to-RGC density'
                    % Set desired cone-RGC density
                    coneToRGCDensityRatio = 5;
        
                    % Instantiation with custom density regular hex RGC mosaic
                    theRGCconnectorOBJ = RGCconnector(theInputConeMosaic, ...
                        'coneToRGCDensityRatio', coneToRGCDensityRatio, ...
                        'chromaticSpatialVarianceTradeoff', chromaticSpatialVarianceTradeoff, ...
                        'rfCentroidOverlapPenaltyFactor', rfCentroidOverlapPenaltyFactor, ...
                        'maxNeighborNormDistance', maxNeighborNormDistance, ...
                        'maxPassesNum', maxPassesNum, ...
                         'maxMeanConeInputsPerRGCToConsiderSwapping', maxMeanConeInputsPerRGCToConsiderSwapping, ...
                        'maxNumberOfConesToSwap', maxNumberOfConesToSwap, ...
                        'visualizeIntermediateConnectivityStages', visualizeIntermediateConnectivityStages);
        
                case 'custom RGC position matrix'
                    % Generate test RGC positions
                    center = theInputConeMosaic.eccentricityMicrons;
                    RFdiameterMicrons = 15;
                    testRGCpositionsMicrons = [];
                    for k = 1:3
                        for l = 1:3
                            kk = size(testRGCpositionsMicrons,1)+1;
                            testRGCpositionsMicrons(kk,1) = center(1) + (k-1)*RFdiameterMicrons;
                            testRGCpositionsMicrons(kk,2) = center(2) + (l-1)*RFdiameterMicrons;
                        end
                    end

                    testRGCspacingsMicrons = zeros(size(testRGCpositionsMicrons,1),1) + RFdiameterMicrons;

                    % Instantiation with custom RGC lattice positions
                    theRGCconnectorOBJ = RGCconnector(theInputConeMosaic, ...
                        'RGCRFpositionsMicrons', testRGCpositionsMicrons, ...
                        'RGCRFspacingsMicrons', testRGCspacingsMicrons, ...
                        'chromaticSpatialVarianceTradeoff', chromaticSpatialVarianceTradeoff, ...
                        'rfCentroidOverlapPenaltyFactor', rfCentroidOverlapPenaltyFactor, ...
                        'maxNeighborNormDistance', maxNeighborNormDistance, ...
                        'maxPassesNum', maxPassesNum, ...
                         'maxMeanConeInputsPerRGCToConsiderSwapping', maxMeanConeInputsPerRGCToConsiderSwapping, ...
                        'maxNumberOfConesToSwap', maxNumberOfConesToSwap, ...
                        'visualizeIntermediateConnectivityStages', visualizeIntermediateConnectivityStages);
            end % switch
        
            toc
        
            % Export final connectivity
            hFig = theRGCconnectorOBJ.visualizeCurrentConnectivityState(9999);
            drawnow
           
            pfdFileName = sprintf('Ecc_%s_MaxNeighborDist_%2.2f_ChromaSpatialVarianceTradeoff_%2.2f_centroidOverlapPenalty_%2.2f.pdf',eccentricity, maxNeighborNormDistance, chromaticSpatialVarianceTradeoff, rfCentroidOverlapPenaltyFactor);
            NicePlot.exportFigToPDF(pfdFileName, hFig, 300);
            
            % Export the generated @RGCconnector object
            theRGCconnectorFileName = strrep(pfdFileName, '.pdf', '.mat');
            save(theRGCconnectorFileName, 'theRGCconnectorOBJ', '-v7.3');
            fprintf('\nGenerated @RGCconnectorOBJ saved in %s\n', theRGCconnectorFileName);
        
        
        end % iSearchIndex
    end % iTradeOff
end

