function testMidgetRGCmosaic

    regenerateMidgetRGCMosaic = true;
    analyzeRetinalRFoverlap = true;
    recomputeRetinalRFoverlap = true;
    reMapRFs = ~true;

    rfOverlapRatio = 0.0;
    horizontalEccsExamined = -[0 1 2 4 6 8 12 16 20 24 30];
    for iEcc = 1:numel(horizontalEccsExamined)
        close all
        eccDegs = [horizontalEccsExamined(iEcc) 0];
        doIt(eccDegs, regenerateMidgetRGCMosaic, analyzeRetinalRFoverlap, recomputeRetinalRFoverlap, reMapRFs, rfOverlapRatio);
    end

end

function doIt(eccDegs, regenerateMidgetRGCMosaic, analyzeRetinalRFoverlap, recomputeRetinalRFoverlap, reMapRFs, rfOverlapRatio)
    
    sizeDegs = [1 1]*max([0.5 0.5+max(abs(eccDegs))])*0.2;

    mosaicFileName = sprintf('theTestMidgetRGCmosaicEcc_%2.1f_%2.1f.mat', eccDegs(1), eccDegs(2));

    if (regenerateMidgetRGCMosaic)
        fprintf('Generating a midgetRGCMosaic with size %2.1f x %2.1f degs\n', sizeDegs(1), sizeDegs(2));
        m = midgetRGCMosaic(...
            'sourceLatticeSizeDegs', 60, ...
            'eccentricityDegs', eccDegs, ...
            'sizeDegs', sizeDegs ...
            );

        save(mosaicFileName, 'm');
    else
        load(mosaicFileName, 'm');
    end

    if (analyzeRetinalRFoverlap)

        hFig = figure(444);
        m.visualizeRFcenterConnectivity('figureHandle', hFig);
        
        nnndFileName = sprintf('NND_%2.2f_Ecc_%2.1f_%2.1f.mat', rfOverlapRatio, eccDegs(1), eccDegs(2));

        if (recomputeRetinalRFoverlap)
            m.adjustRFoverlap(rfOverlapRatio);
        
            % Visualize the retinal RFs
            m.visualizeRetinalRFs('exportGraphicForEachRF', true);
            
            % Analyze the retinal RF overlap
            [NNNDs, NNNDtuplets, RGCdistances, distancesFromMosaicCenterDegs, targetRGCindices] = m.analyzeRetinalRFoverlap();
            
            save(nnndFileName, 'NNNDs', 'NNNDtuplets', 'RGCdistances', 'distancesFromMosaicCenterDegs', 'targetRGCindices');
        else
            load(nnndFileName, 'NNNDs', 'NNNDtuplets', 'RGCdistances', 'distancesFromMosaicCenterDegs', 'targetRGCindices');
        end



        hFig = figure(111);
        set(hFig, 'Color', [1 1 1], 'Position', [300 300 1000 400]);
        clf;
        ax = subplot(1,2,1);
        cMap = [1 0.5 0.5; 0.5 1.0 0.5; 0.5 0.5 1.0];
        cMap2 = [1 0 0; 0 1.0 0; 0 0  1.0];
        for iRGC = 1:size(NNNDtuplets,1)
            ecc = m.rgcRFpositionsDegs(iRGC,1);
            iNeighborIndices = find(~isnan(squeeze(NNNDtuplets(iRGC,:))));
            for ii = 1:numel(iNeighborIndices)
                plot(ax,ecc, NNNDtuplets(iRGC,iNeighborIndices(ii)), ...
                    'o', 'MarkerSize', 10, 'MarkerFaceColor', cMap(ii,:), 'MarkerEdgeColor', cMap2(ii,:));
                hold(ax, 'on');
            end

        end
        set(ax, 'YLim', [0 8], 'YTick', 0:8);
        set(ax, 'FontSize', 16)
        xlabel(ax,'eccentricity (degs)');
        ylabel(ax,'normalized distance');

        ax = subplot(1,2,2);
        edges = 0:0.2:10;
        h = histogram(ax,squeeze(NNNDtuplets(:,1)), edges);
        h.FaceColor = [1 0.5 0.5];
        h.EdgeColor = [1 0 0];
        set(ax, 'FontSize', 16)
        set(ax, 'XLim', [0 8], 'XTick', 0:8);
        title(ax,sprintf('median NNND = %2.2f, ecc: (%2.2f,%2.2f)',median(NNNDs), eccDegs(1), eccDegs(2)));
        xlabel(ax,'normalized nearest neighbor distance');
        ylabel(ax,'cell #');
        NicePlot.exportFigToPDF(strrep(nnndFileName, '.mat', '.pdf'), hFig, 300);
    end


    if (reMapRFs)
        % Generate optics for the mosaic
        mosaicCenterPositionDegs = mean(m.inputConeMosaic.coneRFpositionsDegs,1);
        opticsParams = struct(...
            'positionDegs', mosaicCenterPositionDegs, ...  % (x,y) eccentricity for the PSF, in degrees
            'ZernikeDataBase', 'Artal2012', ...
            'examinedSubjectRankOrder', 2, ...
            'refractiveErrorDiopters', 0.0, ...   % use -999 for optics that do not subtract the central refraction
            'analyzedEye', 'right eye', ...
            'subjectRankingEye', 'right eye', ...
            'pupilDiameterMM', 3.0, ...
            'wavefrontSpatialSamples', 501, ...
            'psfUpsampleFactor', 1 ...
            );
    
        [thePSFData, ~,~, theOI] = RetinaToVisualFieldTransformer.computeVlambdaWeightedPSF(...
            opticsParams, m.inputConeMosaic, m.inputConeMosaic.wave);
    
    
        sceneFOVdegs = m.inputConeMosaic.sizeDegs * 1.1;

        % Generate a presentation display with a desired resolution
        pixelsNum = 256;
        retinalImageResolutionDegs = max(sceneFOVdegs)/pixelsNum;
        viewingDistanceMeters = 4;
        theDisplay = rfMappingStimulusGenerator.presentationDisplay(...
            m.inputConeMosaic.wave, retinalImageResolutionDegs, viewingDistanceMeters);


        % Stim params for the RF mapping
        stimParams = struct(...
            'backgroundLuminanceCdM2', 50.0, ...
            'backgroundChromaticity', [0.301 0.301], ...
            'coneContrasts', [1 1 1], ...
            'contrast', 0.75, ...
            'pixelSizeDegs', retinalImageResolutionDegs, ...
            'stimSizeDegs', max(sceneFOVdegs), ...
            'wavelengthSupport', displayGet(theDisplay, 'wave'), ...
            'viewingDistanceMeters', displayGet(theDisplay, 'viewing distance') ...
            );
    
        % Compute the Hartley spatial patterns
        omega = 11;
        fprintf('Generating Hartley patterns\n')
        % Compute spatial modulation patterns for the Hartley set
        [HartleySpatialModulationPatterns, lIndices, mIndices] = ...
            rfMappingStimulusGenerator.HartleyModulationPatterns(...
            omega, stimParams.stimSizeDegs, stimParams.pixelSizeDegs);


        % Generate scenes for the Hartley patterns
        fprintf('Generating scenes for the Hartley patterns\n')
        [theRFMappingStimulusScenes, theNullStimulusScene, spatialSupportDegs] = ...
            rfMappingStimulusGenerator.generateStimulusMappingFramesOnPresentationDisplay(...
                theDisplay, stimParams, HartleySpatialModulationPatterns, ...
                'validateScenes', false);
    
        % Generate scenes for the inverse-polarity Hartley patterns
        fprintf('Generating scenes for the inverse polarity Hartley patterns\n');
        theInversePolarityRFMappingStimulusScenes = ...
            rfMappingStimulusGenerator.generateStimulusMappingFramesOnPresentationDisplay(...
                theDisplay, stimParams, -HartleySpatialModulationPatterns, ...
                'validateScenes', false);
    
        fprintf('Done with all scene generation\n');
    

        % Preallocate memory for RFs
        rgcsNum = size(m.rgcRFpositionsMicrons,1); 
        nStim = numel(theRFMappingStimulusScenes);
        theRFmaps = zeros(rgcsNum, pixelsNum, pixelsNum);
        theMidgetRGCMosaicResponses = zeros(nStim, rgcsNum);
        theMidgetRGCMosaicInverseResponses = zeros(nStim, rgcsNum);
    
        fprintf('Preallocated memory for RF mapping');

        % Compute the cone mosaic responses and build - up the RF
        for iFrame = 0:nStim
            if (iFrame == 0)
                theOI = oiCompute(theNullStimulusScene, theOI);
                theNullStimMidgetRGCMosaicResponse = m.compute(theOI);
                continue;
            end
    
            fprintf('Computing midget RGC mosaic response to stim %d of %d\n', iFrame, nStim);
            % The forward polarity responses
            theOI = oiCompute(theRFMappingStimulusScenes{iFrame}, theOI);
            theMidgetRGCMosaicResponses(iFrame,:) = ...
                (m.compute(theOI) - theNullStimMidgetRGCMosaicResponse)./theNullStimMidgetRGCMosaicResponse;
    
            % The inverse polarity responses
            theOI = oiCompute(theInversePolarityRFMappingStimulusScenes{iFrame}, theOI);
            theMidgetRGCMosaicInverseResponses(iFrame,:) = ...
                (m.compute(theOI) - theNullStimMidgetRGCMosaicResponse)./theNullStimMidgetRGCMosaicResponse;
    
            % Update the RF map for each RGC
            parfor iRGC = 1:rgcsNum
                r = theMidgetRGCMosaicResponses(iFrame,iRGC) - theMidgetRGCMosaicInverseResponses(iRGC);
                theRFmaps(iRGC,:,:) = theRFmaps(iRGC,:,:) + ...
                    HartleySpatialModulationPatterns(iFrame,:,:) * r;
            end
        end

        theRFmaps = 1/(2*nStim)*theRFmaps;

        save('theTestMidgetRGCmosaic.mat', 'theRFmaps', 'spatialSupportDegs', '-append');
        visualizeRFmaps(m, theRFmaps, spatialSupportDegs);
    end

    
end


function visualizeRFmaps(m, theRFmaps, spatialSupportDegs)
    
    mRGCmosaicCenterDegs = mean(m.rgcRFpositionsDegs,1);
    

    hFig = figure(1);
    set(hFig, 'Color', [1 1 1]);
    clf;

    for iRGC = 1: size(theRFmaps,1)

        [D, idx] = MosaicConnector.pdist2(m.rgcRFpositionsDegs, m.rgcRFpositionsDegs(iRGC,:), ...
            'smallest', 2);
        theClosestNeighborRGC = idx(2);

        
        theRF = squeeze(theRFmaps(iRGC,:,:));
        theRF = theRF / max(theRF(:));
        theNearestRF = squeeze(theRFmaps(theClosestNeighborRGC,:,:));
        theNearestRF = theNearestRF / max(theNearestRF(:));

        theFittedGaussian1 = RetinaToVisualFieldTransformer.fitGaussianEllipsoid(...
            spatialSupportDegs+mRGCmosaicCenterDegs(1), ...
            spatialSupportDegs+mRGCmosaicCenterDegs(1), theRF, ...
            'flatTopGaussian', false, ...
            'forcedOrientationDegs', 0, ...
            'rangeForEllipseRcYRcXratio', [1.0 1.0], ...
            'forcedCentroidXYpos', [], ...
            'globalSearch', true, ...
            'multiStartsNum', 8);

        center1 = theFittedGaussian1.xyCenter;
        sigma1 = theFittedGaussian1.characteristicRadii(1)/sqrt(2.0);
        
        theFittedGaussian2 = RetinaToVisualFieldTransformer.fitGaussianEllipsoid(...
            spatialSupportDegs+mRGCmosaicCenterDegs(1), ...
            spatialSupportDegs+mRGCmosaicCenterDegs(1), theNearestRF, ...
            'flatTopGaussian', false, ...
            'forcedOrientationDegs', 0, ...
            'rangeForEllipseRcYRcXratio', [1.0 1.0], ...
            'forcedCentroidXYpos', [], ...
            'globalSearch', true, ...
            'multiStartsNum', 8);

        center2 = theFittedGaussian2.xyCenter;
        sigma2 = theFittedGaussian2.characteristicRadii(1)/sqrt(2.0);
       
        dd = abs(center1-center2);
        dd2 = dd.^2;
        R = sqrt(sum(dd2(:)));
        NNND(iRGC) = 2*R/(sigma1+sigma2)
    
        ax = subplot(1,3,1);
        x = [-0.2:0.002:0.2];
        g1 = exp(-0.5*(x/sigma1).^2);
        g2 = exp(-0.5*((x-R)/sigma2).^2);
        plot(x, g1, 'r-'); hold on;
        plot(x, g2, 'b-');
        axis 'square';

        ax = subplot(1,3,2);
        % Draw a contour at 1 sigma
        zLevels = exp(-0.5)*[1 1];
        [~,c] = contour(ax,spatialSupportDegs+mRGCmosaicCenterDegs(1), spatialSupportDegs+mRGCmosaicCenterDegs(2), theRF, zLevels);
        c.LineWidth = 1.5;
        c.Color = 'red';
        hold(ax, 'on');
        [~,c] =contour(ax,spatialSupportDegs+mRGCmosaicCenterDegs(1), spatialSupportDegs+mRGCmosaicCenterDegs(2), theNearestRF, zLevels);
        c.LineWidth = 1.5;
        c.Color = 'blue';
        hold(ax, 'off');


    %    connectivityVector = full(squeeze(m.rgcRFcenterConeConnectivityMatrix(:, iRGC)));
    %    inputConeIndices = find(connectivityVector > 0.0001);
    %     for iInput = 1:numel(inputConeIndices)
    %         iCone = inputConeIndices(iInput);
    %         m.inputConeMosaic.coneRFpositionsDegs(iCone,:)
    %         plot(ax, m.inputConeMosaic.coneRFpositionsDegs(iCone,1), m.inputConeMosaic.coneRFpositionsDegs(iCone,2), 'ro');
    %     end
    
    
        set(ax, 'CLim', [0 1]);
        set(ax, 'XLim', [spatialSupportDegs(1) spatialSupportDegs(end)]+mRGCmosaicCenterDegs(1), ...
                'YLim', [spatialSupportDegs(1) spatialSupportDegs(end)]+mRGCmosaicCenterDegs(2));
        axis(ax, 'equal');
        colormap(ax,brewermap(1024, 'greys'));
        
        ax = subplot(1,3,3);
        edges = 1:0.2:4;
        histogram(ax, NNND, edges)
        drawnow
    end


    

end


function compare58Vs60degs()
    computeMosaics = ~true;

    if (~computeMosaics)
        load('midgetRGCs.mat', 'm60_neg25degs', 'm60_pos25degs', 'm58_neg25degs', 'm58_pos25degs');

        m60_neg25degs.inputConeMosaic.horizontalRetinalMeridian
        neg25degsConesNum = numel(m60_neg25degs.inputConeMosaic.coneTypes)
        neg25degsRGCsnum = numel(m60_neg25degs.rgcRFspacingsMicrons)
        
        m60_pos25degs.inputConeMosaic.horizontalRetinalMeridian
        pos25degsConesNum = numel(m60_pos25degs.inputConeMosaic.coneTypes)
        pos25degsRGCsnum = numel(m60_pos25degs.rgcRFspacingsMicrons)

        return;
    end

    m58_pos25degs = midgetRGCMosaic(...
        'sourceLatticeSizeDegs', 58, ...
        'eccentricityDegs', [25 0], ...
        'sizeDegs', [2 2] ...
        );

    m58_neg25degs = midgetRGCMosaic(...
        'sourceLatticeSizeDegs', 58, ...
        'eccentricityDegs', [-25 0], ...
        'sizeDegs', [2 2] ...
        );

    m60_pos25degs = midgetRGCMosaic(...
        'sourceLatticeSizeDegs', 60, ...
        'eccentricityDegs', [25 0], ...
        'sizeDegs', [2 2] ...
        );

    m60_neg25degs = midgetRGCMosaic(...
        'sourceLatticeSizeDegs', 60, ...
        'eccentricityDegs', [-25 0], ...
        'sizeDegs', [2 2] ...
        );

    save('midgetRGCs.mat', 'm60_neg25degs', 'm60_pos25degs', 'm58_neg25degs', 'm58_pos25degs', '-v7.3');
end
