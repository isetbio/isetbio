function measureMRGCRFcentersAcrossXYeccentricityForDifferentSubjects()

    % Horizontal eccentricities examined
    eccX = [-25 -20 -16 -12 -8 -4 0 4 8 12 20 25];
    eccX = [-8];
  

    % Vertical eccentricities examined
    eccY = [0];

    % Optics
    ZernikeDataBase = 'Polans2015';
    subjectRankOrder = 1;

    %newZernikeDataBase = 'Artal2012';
    %newSubjectRankOrder = 3;

    newZernikeDataBase = 'Polans2015';
    newSubjectRankOrder = 10;

    maxVisualizedRFs = 32;

    % Generate the midgetRGCmosaic and optics
    generateTheComponents = ~true;
    if (generateTheComponents)
        generateComponents(ZernikeDataBase, subjectRankOrder, eccX, eccY);
        visualizeComponents(ZernikeDataBase, subjectRankOrder, eccX, eccY, maxVisualizedRFs);
    end

    % Modify the optics
    changeOpticsSubject = true;
    if (changeOpticsSubject)
        modifyOptics(ZernikeDataBase, subjectRankOrder, newZernikeDataBase, newSubjectRankOrder, eccX, eccY);
        visualizeComponents(newZernikeDataBase, newSubjectRankOrder, eccX, eccY, maxVisualizedRFs);
    end
    
    % Compute RF maps (center)
    computeTheVisualRFmaps = true;
    if (computeTheVisualRFmaps)
        computeVisualRF(newZernikeDataBase, newSubjectRankOrder, eccX, eccY, maxVisualizedRFs);
    end

    % Visualize RF maps
    visualizeTheVisualRFmaps = true;
    if (visualizeTheVisualRFmaps)
        visualizeComputedRFs(newZernikeDataBase, newSubjectRankOrder, eccX, eccY);
    end

end

function visualizeComputedRFs(ZernikeDataBase, subjectRankOrder, eccX, eccY)

    [eccXGrid, eccYGrid] = meshgrid(eccX, eccY);
    eccXGrid = eccXGrid(:);
    eccYGrid = eccYGrid(:);

    for iPos = 1:numel(eccXGrid)

        % The position analyzed
        mosaicEccDegs = [eccXGrid(iPos) eccYGrid(iPos)];

        % Load data
        fNameRFmaps = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f_RFmaps.mat', ...
            ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2));

        load(fNameRFmaps, 'theMidgetRGCmosaic', 'opticsParams', 'visualRFspatialSupportDegs', ...
                          'retinalRFcenterMaps', 'theVisualRFmaps', 'rgcIndicesOfAnalyzedRFs');

        [oiEnsemble, psfEnsemble] = theMidgetRGCmosaic.inputConeMosaic.oiEnsembleGenerate(opticsParams.positionDegs, ...
                    'zernikeDataBase', opticsParams.ZernikeDataBase, ...
                    'subjectID', opticsParams.testSubjectID, ...
                    'pupilDiameterMM', opticsParams.pupilDiameterMM, ...
                    'refractiveErrorDiopters', opticsParams.refractiveErrorDiopters, ...
                    'zeroCenterPSF', opticsParams.zeroCenterPSF, ...
                    'subtractCentralRefraction', opticsParams.subtractCentralRefraction, ...
                    'wavefrontSpatialSamples', opticsParams.wavefrontSpatialSamples, ...
                    'upsampleFactor', opticsParams.psfUpsampleFactor, ...
                    'warningInsteadOfErrorForBadZernikeCoeffs', true);
        thePSFData = psfEnsemble{1};
        idx = find(thePSFData.supportWavelength == 550);
        thePSF = thePSFData.data(:,:,idx);

        figure(5);
        ax = subplot(1,3,1);
        imagesc(ax, thePSFData.supportX/60, thePSFData.supportY/60, thePSF/max(thePSF(:)));
        axis(ax,'image'); axis(ax, 'xy');
        set(ax, 'XLim', 0.25*[-1 1], 'YLim', 0.25*[-1 1]);

        xLims = theMidgetRGCmosaic.inputConeMosaic.eccentricityDegs(1) + [-0.25 0.25];
        yLims = theMidgetRGCmosaic.inputConeMosaic.eccentricityDegs(2) + [-0.25 0.25];


        for iRGC = 1:numel(rgcIndicesOfAnalyzedRFs)
            r = retinalRFcenterMaps{iRGC};
            
            ax = subplot(1,3,2);
            imagesc(ax,r.spatialSupportDegsX, r.spatialSupportDegsY, r.centerRF);
            axis(ax,'image'); axis(ax, 'xy');
            set(ax, 'XLim', xLims, 'YLim', yLims);

            ax = subplot(1,3,3);
            imagesc(ax,visualRFspatialSupportDegs+theMidgetRGCmosaic.inputConeMosaic.eccentricityDegs(1), ...
                       visualRFspatialSupportDegs+theMidgetRGCmosaic.inputConeMosaic.eccentricityDegs(2), ...
                       squeeze(theVisualRFmaps(iRGC,:,:)));
            axis(ax,'image'); axis(ax, 'xy');
            set(ax, 'XLim', xLims, 'YLim', yLims, 'Color', [0 0 0]);
            title(ax, sprintf('%s_rank%d', ZernikeDataBase, subjectRankOrder));
            colormap(gray);
            pause(1)
        end

    end % iPos
end


function computeVisualRF(ZernikeDataBase, subjectRankOrder, eccX, eccY, centerMostRGCsNumToAnalyze)
    [eccXGrid, eccYGrid] = meshgrid(eccX, eccY);
    eccXGrid = eccXGrid(:);
    eccYGrid = eccYGrid(:);

    for iPos = 1:numel(eccXGrid)

        % The position analyzed
        mosaicEccDegs = [eccXGrid(iPos) eccYGrid(iPos)];

        % Load the computed components data
        fName = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
            ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2));
        load(fName, 'theMidgetRGCmosaic', 'opticsParams');

        % Compute the selected subject optics using the saved opticsParams
        % Change pupil diameter to 3.5 mm (to match Johannes experiment)
        opticsParams.pupilDiameterMM = 3.0;

        [oiEnsemble, psfEnsemble] = theMidgetRGCmosaic.inputConeMosaic.oiEnsembleGenerate(opticsParams.positionDegs, ...
                    'zernikeDataBase', opticsParams.ZernikeDataBase, ...
                    'subjectID', opticsParams.testSubjectID, ...
                    'pupilDiameterMM', opticsParams.pupilDiameterMM, ...
                    'refractiveErrorDiopters', opticsParams.refractiveErrorDiopters, ...
                    'zeroCenterPSF', opticsParams.zeroCenterPSF, ...
                    'subtractCentralRefraction', opticsParams.subtractCentralRefraction, ...
                    'wavefrontSpatialSamples', opticsParams.wavefrontSpatialSamples, ...
                    'upsampleFactor', opticsParams.psfUpsampleFactor, ...
                    'warningInsteadOfErrorForBadZernikeCoeffs', true);
        theSubjectOptics = oiEnsemble{1};

        % Compute visual RF maps using subspace rev corr
        % Generate a presentation display with a desired resolution
        pixelsNum = 256;
        retinalImageResolutionDegs = 0.5*max(theMidgetRGCmosaic.sizeDegs)/pixelsNum;
        viewingDistanceMeters = 4;
        theDisplay = rfMappingStimulusGenerator.presentationDisplay(...
            theMidgetRGCmosaic.inputConeMosaic.wave, retinalImageResolutionDegs, ...
            viewingDistanceMeters);

        % Stim params for the RF mapping
        stimParams = struct(...
            'backgroundLuminanceCdM2', 50.0, ...
            'backgroundChromaticity', [0.301 0.301], ...
            'coneContrasts', [1 1 1], ...
            'contrast', 0.75, ...
            'pixelSizeDegs', retinalImageResolutionDegs, ...
            'stimSizeDegs', 0.5*max(theMidgetRGCmosaic.sizeDegs), ...  % make the stimulus size = 1/2 * RGC mosaic FoV
            'wavelengthSupport', displayGet(theDisplay, 'wave'), ...
            'viewingDistanceMeters', displayGet(theDisplay, 'viewing distance') ...
            );

        % Hartley (RF mapping) spatial patterns
        fprintf('Generating Hartley patterns\n');
        omega = 9;
        % Compute spatial modulation patterns for the Hartley set
        HartleySpatialModulationPatterns = ...
            rfMappingStimulusGenerator.HartleyModulationPatterns(...
            omega, stimParams.stimSizeDegs, stimParams.pixelSizeDegs);
    

        
        % Find the indices of the centerMostRGCsNumToAnalyze RGCs
        relativeRGCpositions = bsxfun(@minus, theMidgetRGCmosaic.rgcRFpositionsDegs, theMidgetRGCmosaic.inputConeMosaic.eccentricityDegs);
        radii = sum(relativeRGCpositions.^2,2);
        [~,sortedRGCindices] = sort(radii, 'ascend');
        rgcsNum = size(theMidgetRGCmosaic.rgcRFpositionsDegs,1);
        rgcIndicesOfAnalyzedRFs = sortedRGCindices(1:min([centerMostRGCsNumToAnalyze, rgcsNum]));

        % Preallocate memory
        nStim = size(HartleySpatialModulationPatterns,1);
        pixelsNum = size(HartleySpatialModulationPatterns,2);
        theMRGCMosaicResponses = zeros(nStim, numel(rgcIndicesOfAnalyzedRFs));
        theMRGCMosaicReversePolarityResponse = zeros(nStim, numel(rgcIndicesOfAnalyzedRFs));

        % Compute theNullStimulusScene and the visual RF spatial support
        [~, theNullStimulusScene, visualRFspatialSupportDegs] = ...
            rfMappingStimulusGenerator.generateStimulusMappingFramesOnPresentationDisplay(...
                    theDisplay, stimParams, HartleySpatialModulationPatterns(1,:,:), ...
                    'validateScenes', false);

        % Compute the MRGCmosaic responses to the forward polarity stimuli
        fprintf('Computing MRGC responses to the positive polarity stimuli ...\n');

        parfor iFrame = 1:nStim
            fprintf('Forward polarity stimulus %d/%d\n', iFrame, nStim);

            % Generate scene for the Hartley pattern
            theRFMappingStimulusFrameScene = ...
                rfMappingStimulusGenerator.generateStimulusMappingFramesOnPresentationDisplay(...
                    theDisplay, stimParams, HartleySpatialModulationPatterns(iFrame,:,:), ...
                    'validateScenes', false);

            % Compute the mosaic's response to the positive polarity frame
            r = theMidgetRGCmosaic.compute(...
                    theRFMappingStimulusFrameScene{1}, ...
                    'withOptics', theSubjectOptics, ...
                    'nTrials', 1, ...
                    'theNullScene', theNullStimulusScene, ...
                    'normalizeConeResponsesWithRespectToNullScene', true);

            % Save response to this stimulus
            theMRGCMosaicResponses(iFrame,:) = (squeeze(r(rgcIndicesOfAnalyzedRFs)))';
        end


        % Compute the MRGCmosaic responses to the reverse polarity stimuli
        fprintf('Computing MRGC responses to the reverse polarity stimuli ...\n');

        parfor iFrame = 1:nStim
            fprintf('Reverse polarity stimulus %d/%d\n', iFrame, nStim);

            % Generate scene for the Hartley pattern
            theRFMappingStimulusFrameScene = ...
                rfMappingStimulusGenerator.generateStimulusMappingFramesOnPresentationDisplay(...
                    theDisplay, stimParams, -HartleySpatialModulationPatterns(iFrame,:,:), ...
                    'validateScenes', false);

            % Compute the mosaic's response to the reverse polarity frames
            r = theMidgetRGCmosaic.compute(...
                    theRFMappingStimulusFrameScene{1}, ...
                    'withOptics', theSubjectOptics, ...
                    'nTrials', 1, ...
                    'theNullScene', theNullStimulusScene, ...
                    'normalizeConeResponsesWithRespectToNullScene', true);

             % Save response to this stimulus
             theMRGCMosaicReversePolarityResponse(iFrame,:) = (squeeze(r(rgcIndicesOfAnalyzedRFs)))';
        end

        % Compute forward - reverse polarity responses
        theMRGCMosaicResponses = theMRGCMosaicResponses - theMRGCMosaicReversePolarityResponse;

        clear 'theMRGCMosaicReversePolarityResponse'

        % Allocate memory for the visual RF maps (single to save space)
        theVisualRFmaps = zeros(numel(rgcIndicesOfAnalyzedRFs), pixelsNum, pixelsNum, 'single');

        % Also make the Hartley patterns and the responses single
        theMRGCMosaicResponses = single(theMRGCMosaicResponses);

        % Compute the RFs
        parfor iRGCindex = 1:numel(rgcIndicesOfAnalyzedRFs)
            for iFrame = 1:nStim
                theVisualRFmaps(iRGCindex,:,:) = theVisualRFmaps(iRGCindex,:,:) + ...
                    HartleySpatialModulationPatterns(iFrame,:,:) * theMRGCMosaicResponses(iFrame,iRGCindex);
            end
        end % iRGCindex

        % Normalize
        theVisualRFmaps = 1/(2*nStim)*theVisualRFmaps;

        % Compute retinal RFmaps
        marginDegs = min([0.5 0.4*min(theMidgetRGCmosaic.sizeDegs)]);
        spatialSupportSamplesNum = 256;
        retinalRFcenterMaps = theMidgetRGCmosaic.computeRetinalRFcenterMaps(...
            marginDegs, spatialSupportSamplesNum, ...
            'forRGCindices', rgcIndicesOfAnalyzedRFs);

        for iRGCindex = 1:numel(retinalRFcenterMaps)
            r = retinalRFcenterMaps{iRGCindex};
            retinalRFcenterMap = r.centerRF;
            retinalSpatialSupportDegsX = r.spatialSupportDegsX;
            retinalSpatialSupportDegsY = r.spatialSupportDegsY;

            xLims = [min(retinalSpatialSupportDegsX ) max(retinalSpatialSupportDegsX )];
            yLims = [min(retinalSpatialSupportDegsY) max(retinalSpatialSupportDegsY)];

            cMap = gray(1024);
            hFig = figure(1); clf;
            set(hFig, 'Position', [10 10 1700 850], 'Color', [1 1 1]);
            subplot(1,2,1);
            imagesc(retinalSpatialSupportDegsX, retinalSpatialSupportDegsY, retinalRFcenterMap);
            axis 'image'; axis 'xy';
            set(gca, 'XLim', xLims, 'YLim', yLims, 'FontSize', 20);
            title(sprintf('Retinal RF map \nX,Y = (%2.1f,%2.1f) degs', ...
                theMidgetRGCmosaic.rgcRFpositionsDegs(rgcIndicesOfAnalyzedRFs(iRGCindex),1), theMidgetRGCmosaic.rgcRFpositionsDegs(rgcIndicesOfAnalyzedRFs(iRGCindex),2)));
            
            subplot(1,2,2);
            theRFmap = squeeze(theVisualRFmaps(iRGCindex,:,:));
            theRFmap = theRFmap / max(abs(theRFmap(:)));

            imagesc(visualRFspatialSupportDegs+theMidgetRGCmosaic.inputConeMosaic.eccentricityDegs(1), ...
                    visualRFspatialSupportDegs+theMidgetRGCmosaic.inputConeMosaic.eccentricityDegs(2), ...
                    theRFmap);
            axis 'image'; axis 'xy';
            set(gca, 'XLim', xLims, 'YLim', yLims, 'FontSize', 20, 'Color', cMap(512,:), 'CLim', [-1 1]);
            title(sprintf('Visual RF map \nX,Y = (%2.1f,%2.1f) degs', ...
                theMidgetRGCmosaic.rgcRFpositionsDegs(rgcIndicesOfAnalyzedRFs(iRGCindex),1), theMidgetRGCmosaic.rgcRFpositionsDegs(rgcIndicesOfAnalyzedRFs(iRGCindex),2)));
            colormap(cMap);
        end % iRGCindex

        % Save data
        fNameRFmaps = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f_RFmaps.mat', ...
            ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2));

        save(fNameRFmaps, 'theMidgetRGCmosaic', 'opticsParams', 'visualRFspatialSupportDegs', ...
                          'retinalRFcenterMaps', 'theVisualRFmaps', 'rgcIndicesOfAnalyzedRFs', '-v7.3');
        fprintf('RFmaps for subject %d at position %d of %d saved to %s\n', ...
            subjectRankOrder, iPos, numel(eccXGrid), fNameRFmaps);

    end % iPos

end

function  visualizeComponents(ZernikeDataBase, subjectRankOrder, eccX, eccY, maxVisualizedRFs)

    [eccXGrid, eccYGrid] = meshgrid(eccX, eccY);
    eccXGrid = eccXGrid(:);
    eccYGrid = eccYGrid(:);

    rowsNum = 1; %numel(eccX);
    colsNum = 1; %numel(eccY);
    sv = NicePlot.getSubPlotPosVectors(...
       'colsNum', colsNum, ...
       'rowsNum', rowsNum, ...
       'heightMargin',  0.04, ...
       'widthMargin',    0.01, ...
       'leftMargin',     0.10, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.08, ...
       'topMargin',      0.02); 

    for iPos = 1:numel(eccXGrid)

        hFig = figure(1);
        set(hFig, 'Position', [10 10 900 900], 'Color', [1 1 1]);

        % The position analyzed
        mosaicEccDegs = [eccXGrid(iPos) eccYGrid(iPos)];

        % Load the computed components data
        fName = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
            ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2));

        load(fName, 'thePSFData', 'theMidgetRGCmosaic', 'opticsParams');
        
        r = floor((iPos-1)/colsNum);
        r = mod(r,rowsNum)+1;
        c = mod(iPos-1,colsNum)+1;
        ax = subplot('Position', sv(r,c).v);

        theMidgetRGCmosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax, ...
            'maxVisualizedRFs', maxVisualizedRFs, ...
            'xRange', 0.75, ...
            'yRange', 0.75, ...
            'fontSize', 20);

        % Overlay the PSF
        hold(ax, 'on');
        cmap = brewermap(1024,'greys');
        alpha = 0.75;
        contourLineColor = [0 0 0.0];
        mRGCmosaicCenterDegs = mean(theMidgetRGCmosaic.rgcRFpositionsDegs,1);
        cMosaic.semiTransparentContourPlot(ax, ...
            thePSFData.psfSupportXdegs + mRGCmosaicCenterDegs(1), ...
            thePSFData.psfSupportYdegs + mRGCmosaicCenterDegs(2), ...
            thePSFData.vLambdaWeightedData/max(thePSFData.vLambdaWeightedData(:)), ...
            0.05:0.15:0.95, cmap, alpha, contourLineColor, ...
            'lineWidth', 2.0);

        fNamePDF = strrep(fName, '.mat', '.pdf');
        NicePlot.exportFigToPDF(fNamePDF, hFig, 300);
    end

end

function modifyOptics(ZernikeDataBase, subjectRankOrder, newZernikeDataBase, newSubjectRankOrder, eccX, eccY)
    % Struct with the various optics params
    opticsParams = struct(...
        'positionDegs', [], ...           % (x,y) eccentricity for the PSF, in degrees
        'ZernikeDataBase', ZernikeDataBase, ...
        'examinedSubjectRankOrder', [], ...
        'refractiveErrorDiopters', 0.0, ...    % use -999 for optics that do not subtract the central refraction
        'analyzedEye', 'right eye', ...
        'subjectRankingEye', 'right eye', ...
        'pupilDiameterMM', 3.0, ...
        'wavefrontSpatialSamples', 701, ...
        'psfUpsampleFactor', 1 ...
        );


    [eccXGrid, eccYGrid] = meshgrid(eccX, eccY);
    eccXGrid = eccXGrid(:);
    eccYGrid = eccYGrid(:);

    for iPos = 1:numel(eccXGrid)

        % The position analyzed
        mosaicEccDegs = [eccXGrid(iPos) eccYGrid(iPos)];
        radialEccDegs = sqrt(sum(mosaicEccDegs.^2,2));

        if (radialEccDegs == 0)
            mosaicSizeDegs = 0.3;
        elseif (radialEccDegs <= 2)
            mosaicSizeDegs = 0.4;
        elseif (radialEccDegs <= 3)
            mosaicSizeDegs = 0.6;
        elseif (radialEccDegs <= 4)
            mosaicSizeDegs = 0.8;
        elseif (radialEccDegs <= 6)
            mosaicSizeDegs = 1.0;
        elseif (radialEccDegs <= 8)
            mosaicSizeDegs = 1.4;
        elseif (radialEccDegs <= 10)
            mosaicSizeDegs = 1.5;
        elseif (radialEccDegs <= 12)
            mosaicSizeDegs = 1.6;
        elseif (radialEccDegs <= 14)
            mosaicSizeDegs = 1.8;
        elseif (radialEccDegs <= 16)
            mosaicSizeDegs = 2.0;
        elseif (radialEccDegs <= 20)
            mosaicSizeDegs = 2.2;
        else
            mosaicSizeDegs = 2.5;
        end

        % Change position
        opticsParams.positionDegs = mosaicEccDegs;

        fName = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
            ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2));
        load(fName, 'theMidgetRGCmosaic');


        % Generate PSF for new subject
        opticsParams.examinedSubjectRankOrder = newSubjectRankOrder;
        opticsParams.ZernikeDataBase = newZernikeDataBase;

        [thePSFData,~,~,opticsParams] = RetinaToVisualFieldTransformer.computeVlambdaWeightedPSF(...
            opticsParams, theMidgetRGCmosaic.inputConeMosaic, []);


        % Save data
        fName = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
            newZernikeDataBase, newSubjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2));

        save(fName, 'thePSFData', 'theMidgetRGCmosaic', 'opticsParams', '-v7.3');
        fprintf('Data for subject %d at position %d of %d saved to %s\n', newSubjectRankOrder, iPos, numel(eccXGrid), fName);
    end

end


function generateComponents(ZernikeDataBase, subjectRankOrder, eccX, eccY)
    % Struct with the various optics params
    opticsParams = struct(...
        'positionDegs', [], ...           % (x,y) eccentricity for the PSF, in degrees
        'ZernikeDataBase', ZernikeDataBase, ...
        'examinedSubjectRankOrder', subjectRankOrder, ...
        'refractiveErrorDiopters', 0.0, ...    % use -999 for optics that do not subtract the central refraction
        'analyzedEye', 'right eye', ...
        'subjectRankingEye', 'right eye', ...
        'pupilDiameterMM', 3.0, ...
        'wavefrontSpatialSamples', 701, ...
        'psfUpsampleFactor', 1 ...
        );


    [eccXGrid, eccYGrid] = meshgrid(eccX, eccY);
    eccXGrid = eccXGrid(:);
    eccYGrid = eccYGrid(:);

    for iPos = 1:numel(eccXGrid)

        % The position analyzed
        mosaicEccDegs = [eccXGrid(iPos) eccYGrid(iPos)];
        radialEccDegs = sqrt(sum(mosaicEccDegs.^2,2));

        if (radialEccDegs == 0)
            mosaicSizeDegs = 0.3;
        elseif (radialEccDegs <= 2)
            mosaicSizeDegs = 0.4;
        elseif (radialEccDegs <= 3)
            mosaicSizeDegs = 0.6;
        elseif (radialEccDegs <= 4)
            mosaicSizeDegs = 0.8;
        elseif (radialEccDegs <= 6)
            mosaicSizeDegs = 1.0;
        elseif (radialEccDegs <= 8)
            mosaicSizeDegs = 1.4;
        elseif (radialEccDegs <= 10)
            mosaicSizeDegs = 1.5;
        elseif (radialEccDegs <= 12)
            mosaicSizeDegs = 1.6;
        elseif (radialEccDegs <= 14)
            mosaicSizeDegs = 1.8;
        elseif (radialEccDegs <= 16)
            mosaicSizeDegs = 2.0;
        elseif (radialEccDegs <= 20)
            mosaicSizeDegs = 2.2;
        else
            mosaicSizeDegs = 2.5;
        end

        % Change position
        opticsParams.positionDegs = mosaicEccDegs;

        % Generate mRGC mosaic
        theMidgetRGCmosaic = midgetRGCMosaic(...
                        'sourceLatticeSizeDegs', 60, ...
                        'whichEye', opticsParams.analyzedEye, ...
                        'eccentricityDegs', mosaicEccDegs, ...
                        'sizeDegs', mosaicSizeDegs*[1 1] ...
                        );

        % Generate PSF
        [thePSFData, ~, ~, opticsParams] = RetinaToVisualFieldTransformer.computeVlambdaWeightedPSF(...
            opticsParams, theMidgetRGCmosaic.inputConeMosaic, []);

        % Save data
        fName = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
            ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2));

        save(fName, 'thePSFData', 'theMidgetRGCmosaic', 'opticsParams', '-v7.3');
        fprintf('Data for position %d of %d saved to %s\n', iPos, numel(eccXGrid), fName);
    end

end
