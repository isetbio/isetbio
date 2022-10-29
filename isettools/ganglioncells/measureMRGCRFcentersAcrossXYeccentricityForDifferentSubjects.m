function measureMRGCRFcentersAcrossXYeccentricityForDifferentSubjects()

    % Horizontal eccentricities examined
    eccX = [-25 -20 -16 -12 -8 -4 0 4 8 12 20 25];
    eccX = [-8];
  

    % Vertical eccentricities examined
    eccX = [-25 -20 -16 -12 -10 -8:1:8 10 12 16 20 25];
    eccY = [-8:1:8];

    eccX = [0];
    eccY = [0];

    % Optics
    ZernikeDataBase = 'Polans2015';
    subjectRankOrder = 1;

    %newZernikeDataBase = 'Artal2012';
    %newSubjectRankOrder = 3;

    newZernikeDataBase = 'Polans2015';
    newSubjectRankOrder = 6;
    newPupilDiamMM = 3.5;

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
        modifyOptics(ZernikeDataBase, subjectRankOrder, newZernikeDataBase, newSubjectRankOrder, newPupilDiamMM, eccX, eccY);
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

    opticalImagePositionDegs = 'mosaic-centered';  % [0.05 -0.03];  %'mosaic-centered'

    visualizedOpticalImageNum = 1;


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
        retinalImageResolutionDegs = max(theMidgetRGCmosaic.sizeDegs)/pixelsNum;
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
            'stimSizeDegs', max(theMidgetRGCmosaic.sizeDegs), ...  % make the stimulus size = 1/2 * RGC mosaic FoV
            'wavelengthSupport', displayGet(theDisplay, 'wave'), ...
            'viewingDistanceMeters', displayGet(theDisplay, 'viewing distance') ...
            );

        % Hartley (RF mapping) spatial patterns
        fprintf('Generating Hartley patterns\n');
        omega = 11;
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
        theVisualRFmaps = zeros(numel(rgcIndicesOfAnalyzedRFs), pixelsNum, pixelsNum, 'single');

        % Compute theNullStimulusScene and the visual RF spatial support
        [~, theNullStimulusScene, visualRFspatialSupportDegs] = ...
            rfMappingStimulusGenerator.generateStimulusMappingFramesOnPresentationDisplay(...
                    theDisplay, stimParams, HartleySpatialModulationPatterns(1,:,:), ...
                    'validateScenes', false);


        % Compute the MRGCmosaic responses to the forward polarity stimuli
        fprintf('Computing MRGC responses to the positive polarity stimuli ...\n');

        for iFrame = 1:nStim
            fprintf('Computing responses and RFs for stimulus %d/%d\n', iFrame, nStim);

            theHartleyPattern = HartleySpatialModulationPatterns(iFrame,:,:);

            % Generate scene for the forward Hartley pattern
            theRFMappingStimulusFrameScene = ...
                rfMappingStimulusGenerator.generateStimulusMappingFramesOnPresentationDisplay(...
                    theDisplay, stimParams, theHartleyPattern, ...
                    'validateScenes', false);

            % Compute the mosaic's response to the positive polarity frame
            [rPlus, ~, noiseFreeConeAbsorptionsCount, theOpticalImage] = theMidgetRGCmosaic.compute(...
                    theRFMappingStimulusFrameScene{1}, ...
                    'withOptics', theSubjectOptics, ...
                    'nTrials', 1, ...
                    'opticalImagePositionDegs', opticalImagePositionDegs, ...
                    'theNullScene', theNullStimulusScene, ...
                    'normalizeConeResponsesWithRespectToNullScene', true);

            % Generate scene for the inverse Hartley pattern
            theRFMappingStimulusFrameScene = ...
                rfMappingStimulusGenerator.generateStimulusMappingFramesOnPresentationDisplay(...
                    theDisplay, stimParams, -theHartleyPattern, ...
                    'validateScenes', false);

            % Compute the mosaic's response to the reverse polarity frames
            rMinus = theMidgetRGCmosaic.compute(...
                    theRFMappingStimulusFrameScene{1}, ...
                    'withOptics', theSubjectOptics, ...
                    'nTrials', 1, ...
                    'opticalImagePositionDegs', opticalImagePositionDegs, ...
                    'theNullScene', theNullStimulusScene, ...
                    'normalizeConeResponsesWithRespectToNullScene', true);

            % Save differentual response to this stimulus frame
            diffResponse = rPlus-rMinus;
            diffResponse = squeeze(diffResponse(rgcIndicesOfAnalyzedRFs));

            % Accumulate the RFs
            parfor iRGCindex = 1:numel(rgcIndicesOfAnalyzedRFs)
                theVisualRFmaps(iRGCindex,:,:) = theVisualRFmaps(iRGCindex,:,:) + ...
                    theHartleyPattern * diffResponse(iRGCindex);
            end

            % Visualize first frame
            if (~isempty(visualizedOpticalImageNum)) && (visualizedOpticalImageNum  == iFrame)
                
                hFig = figure(75);

                % Visualize the input cone mosaic with the OI
                ax = subplot(2,2,1);
                theMidgetRGCmosaic.inputConeMosaic.visualize(...
                    'figureHandle', hFig, ...
                    'axesHandle', ax, ...
                    'withSuperimposedOpticalImage', theOpticalImage);

                % Visualize the input cone mosaic response to the OI
                ax = subplot(2,2,2);
                theMidgetRGCmosaic.inputConeMosaic.visualize(...
                    'figureHandle', hFig, ...
                    'axesHandle', ax, ...
                    'activation', noiseFreeConeAbsorptionsCount);

                % Visualize the midgetRCCmosaic with the OI
                ax = subplot(2,2,3);
                theMidgetRGCmosaic.visualize(...
                    'figureHandle', hFig, ...
                    'axesHandle', ax, ...
                    'withSuperimposedOpticalImage', theOpticalImage);
                drawnow;
                pause
            end

            
        end  % iFrame


        % Normalize the RFs
        theVisualRFmaps = 1/(2*nStim)*theVisualRFmaps;

        % Compute retinal RFmaps
        marginDegs = min([0.5 0.4*min(theMidgetRGCmosaic.sizeDegs)]);
        spatialSupportSamplesNum = 256;
        retinalRFcenterMaps = theMidgetRGCmosaic.computeRetinalRFcenterMaps(...
            marginDegs, spatialSupportSamplesNum, ...
            'forRGCindices', rgcIndicesOfAnalyzedRFs);

        xLims = theMidgetRGCmosaic.eccentricityDegs(1) + theMidgetRGCmosaic.sizeDegs(1)*0.5*[-1 1];
        yLims = theMidgetRGCmosaic.eccentricityDegs(2) + theMidgetRGCmosaic.sizeDegs(2)*0.5*[-1 1];

        for iRGCindex = 1:numel(retinalRFcenterMaps)
            r = retinalRFcenterMaps{iRGCindex};
            retinalRFcenterMap = r.centerRF;
            retinalSpatialSupportDegsX = r.spatialSupportDegsX;
            retinalSpatialSupportDegsY = r.spatialSupportDegsY;

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

        psfData.supportXdegs = thePSFData.psfSupportXdegs;
        psfData.supportYdegs = thePSFData.psfSupportYdegs;
        psfData.data = thePSFData.vLambdaWeightedData;

        theMidgetRGCmosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax, ...
            'maxVisualizedRFs', maxVisualizedRFs, ...
            'withSuperimposedPSF', psfData, ...
            'xRange', 0.75, ...
            'yRange', 0.75, ...
            'fontSize', 20);

        fNamePDF = strrep(fName, '.mat', '.pdf');
        NicePlot.exportFigToPDF(fNamePDF, hFig, 300);
    end

end

function modifyOptics(ZernikeDataBase, subjectRankOrder, newZernikeDataBase, newSubjectRankOrder, newPupilDiamMM, eccX, eccY)
    % Struct with the various optics params
    opticsParams = struct(...
        'positionDegs', [], ...           % (x,y) eccentricity for the PSF, in degrees
        'ZernikeDataBase', ZernikeDataBase, ...
        'examinedSubjectRankOrder', [], ...
        'refractiveErrorDiopters', 0.0, ...    % use -999 for optics that do not subtract the central refraction
        'analyzedEye', 'right eye', ...
        'subjectRankingEye', 'right eye', ...
        'pupilDiameterMM', newPupilDiamMM, ...
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
