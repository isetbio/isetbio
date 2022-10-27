function measureMRGCRFcentersAcrossXYeccentricityForDifferentSubjects()

    % Horizontal eccentricities examined
    eccX = [-24 -20 -16 -12 -8 -4 0 4 8 12 16 20 24];
    eccX = [-25 25];

    % Vertical eccentricities examined
    eccY = [0];

    % Optics
    ZernikeDataBase = 'Polans2015';
    subjectRankOrder = 3;

    generateTheComponents = ~true;
    % Generate the components
    if (generateTheComponents)
        generateComponents(ZernikeDataBase, subjectRankOrder, eccX, eccY);
        visualizeComponents(ZernikeDataBase, subjectRankOrder, eccX, eccY);
    end

    changeOpticsSubject = ~true;
    if (changeOpticsSubject)
        newSubjectRankOrder = 3;
        modifyOptics(ZernikeDataBase, newSubjectRankOrder, eccX, eccY);
        visualizeComponents(ZernikeDataBase, newSubjectRankOrder, eccX, eccY);
    end

    computeVisualRFmap = true;
    if (computeVisualRFmap)
        computeVisualRF(ZernikeDataBase, subjectRankOrder, eccX, eccY);
    end


end

function computeVisualRF(ZernikeDataBase, subjectRankOrder, eccX, eccY)
    [eccXGrid, eccYGrid] = meshgrid(eccX, eccY);
    eccXGrid = eccXGrid(:);
    eccYGrid = eccYGrid(:);

    for iPos = 1:numel(eccXGrid)

        % The position analyzed
        mosaicEccDegs = [eccXGrid(iPos) eccYGrid(iPos)];

        % Load the computed components data
        fName = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
            ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2));
        load(fName, 'thePSFData', 'theMidgetRGCmosaic', 'opticsParams');

        % Compute the OI
        opticsParams
        oiEnsemble = theMidgetRGCmosaic.inputConeMosaic.oiEnsembleGenerate(opticsParams.positionDegs, ...
                   'zernikeDataBase', opticsParams.ZernikeDataBase, ...
                    'subjectID', opticsParams.examinedSubjectRankOrder, ...
                    'pupilDiameterMM', opticsParams.pupilDiameterMM, ...
                    'refractiveErrorDiopters', opticsParams.refractiveErrorDiopters, ...
                    'zeroCenterPSF', opticsParams.zeroCenterPSF, ...
                    'subtractCentralRefraction', opticsParams.subtractCentralRefraction, ...
                    'wavefrontSpatialSamples', opticsParams.wavefrontSpatialSamples, ...
                    'upsampleFactor', opticsParams.psfUpsampleFactor, ...
                    'warningInsteadOfErrorForBadZernikeCoeffs', true);
        theSelectedSubjectOptics = oiEnsemble{1};

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
            'stimSizeDegs', max(theMidgetRGCmosaic.inputConeMosaic.sizeDegs), ...  % make the stimulus = coneMosaicSize
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
    
        % Generate scenes for the Hartley patterns
        fprintf('Generating scenes for the positive polarity stimuli ...\n');
        [theRFMappingStimulusScenes, theNullStimulusScene, visualRFspatialSupportDegs] = rfMappingStimulusGenerator.generateStimulusMappingFramesOnPresentationDisplay(...
            theDisplay, stimParams, HartleySpatialModulationPatterns, ...
            'validateScenes', false);

        % Preallocate memory
        rgcsNum = size(theMidgetRGCmosaic.rgcRFpositionsDegs,1);
        nStim = numel(theRFMappingStimulusScenes);
        pixelsNum = size(HartleySpatialModulationPatterns,2);
        theMRGCMosaicResponses = zeros(nStim, rgcsNum);
        theVisualRFmaps = zeros(rgcsNum, pixelsNum, pixelsNum);

        fprintf('Computing MRGC responses to the positive polarity stimuli ...\n')
        % Compute the MRGCmosaic responses and build - up the RF
        for iFrame = 1:nStim
            % Compute the mosaic's response to the positive polarity frame
            theMRGCMosaicResponses(iFrame,:) = theMidgetRGCmosaic.compute(...
                    theRFMappingStimulusScenes{iFrame}, ...
                    'withOptics', theSelectedSubjectOptics, ...
                    'nTrials', 1, ...
                    'theNullScene', theNullStimulusScene, ...
                    'normalizeConeResponsesWithRespectToNullScene', true);
            % Release some memory
            theRFMappingStimulusScenes{iFrame} = [];
        end

        % Release memory
        clear 'theRFMappingStimulusScenes';

        fprintf('Generating scenes for the reverse polarity stimuli ... \n')
        % Generate scenes for the inverse-polarity Hartley patterns
        theReversePolarityRFMappingStimulusScenes = rfMappingStimulusGenerator.generateStimulusMappingFramesOnPresentationDisplay(...
            theDisplay, stimParams, -HartleySpatialModulationPatterns, ...
            'validateScenes', false);

        fprintf('Computing MRGC responses to the reverse polarity stimuli ...\n');
        for iFrame = 1:nStim
            % Compute the mosaic's response to the reverse polarity frames
            theMRGCMosaicReversePolarityResponse = theMidgetRGCmosaic.compute(...
                    theReversePolarityRFMappingStimulusScenes{iFrame}, ...
                    'withOptics', theSelectedSubjectOptics, ...
                    'nTrials', 1, ...
                    'theNullScene', theNullStimulusScene, ...
                    'normalizeConeResponsesWithRespectToNullScene', true);

            % Release some memory
            theReversePolarityRFMappingStimulusScenes{iFrame} = [];

            % Update the RF map for each cone
            parfor iRGCindex = 1:rgcsNum 
                r = theMRGCMosaicResponses(iFrame,iRGCindex) - theMRGCMosaicReversePolarityResponse(iRGCindex);
                theVisualRFmaps(iRGCindex,:,:) = theVisualRFmaps(iRGCindex,:,:) + ...
                    HartleySpatialModulationPatterns(iFrame,:,:) * r;
            end

        end % iFrame

        % Release memory
        clear 'theReversePolarityRFMappingStimulusScenes'

        % Normalize
        theVisualRFmaps = 1/(2*nStim)*theVisualRFmaps;


        % Compute retinal RFmaps
        marginDegs = min([0.5 0.4*min(theMidgetRGCmosaic.sizeDegs)]);
        spatialSupportSamplesNum = 256;
        retinalRFcenterMaps = theMidgetRGCmosaic.computeRetinalRFcenterMaps(...
            marginDegs, spatialSupportSamplesNum, ...
            'forRGCindices', 1:10);

        for iRGCindex = 1:numel(retinalRFcenterMaps)
            r = retinalRFcenterMaps{iRGCindex};
            retinalRFcenterMap = r.centerRF;
            retinalSpatialSupportDegsX = r.spatialSupportDegsX;
            retinalSpatialSupportDegsY = r.spatialSupportDegsY;
            figure(1); clf;
            subplot(1,2,1);
            imagesc(retinalSpatialSupportDegsX, retinalSpatialSupportDegsY, retinalRFcenterMap);
            axis 'image'
            title(sprintf('Retinal RF map \nX,Y = (%2.1f,%2.1f) degs', ...
                theMidgetRGCmosaic.rgcRFpositionsDegs(1), theMidgetRGCmosaic.rgcRFpositionsDegs(2)));

            subplot(1,2,2);
            imagesc(visualRFspatialSupportDegs, visualRFspatialSupportDegs, squeeze(theVisualRFmaps(iRGCindex,:,:)));
            axis 'image'
            title(sprintf('Visual RF map \nX,Y = (%2.1f,%2.1f) degs', ...
                theMidgetRGCmosaic.rgcRFpositionsDegs(1), theMidgetRGCmosaic.rgcRFpositionsDegs(2)));

        end % iRGCindex

        % Save data
        fNameRFmaps = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f_RFmaps.mat', ...
            ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2));

        save(fNameRFmaps, 'retinalSpatialSupportDegsX', 'retinalSpatialSupportDegsY', 'visualRFspatialSupportDegs', ...
            'retinalRFcenterMaps', 'theVisualRFmaps', '-v7.3');
        fprintf('RFmaps for subject %d at position %d of %d saved to %s\n', subjectRankOrder, iPos, numel(eccXGrid), fNameRFmaps);

    end % iPos

end

function  visualizeComponents(ZernikeDataBase, subjectRankOrder, eccX, eccY)

    [eccXGrid, eccYGrid] = meshgrid(eccX, eccY);
    eccXGrid = eccXGrid(:)
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
            'maxVisualizedRFs', 32, ...
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
            thePSFData.data/max(thePSFData.data(:)), ...
            0.05:0.15:0.95, cmap, alpha, contourLineColor, ...
            'lineWidth', 2.0);

        fNamePDF = strrep(fName, '.mat', '.pdf');
        NicePlot.exportFigToPDF(fNamePDF, hFig, 300);
    end

end

function modifyOptics(ZernikeDataBase, newSubjectRankOrder, eccX, eccY)
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
            mosaicSizeDegs = 1.2;
        elseif (radialEccDegs <= 10)
            mosaicSizeDegs = 1.4;
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
            ZernikeDataBase, 3, mosaicEccDegs(1), mosaicEccDegs(2));
        load(fName, 'theMidgetRGCmosaic');


        % Generate PSF for new subject
        opticsParams.examinedSubjectRankOrder = newSubjectRankOrder;
        [thePSFData,~,~,opticsParams] = RetinaToVisualFieldTransformer.computeVlambdaWeightedPSF(...
            opticsParams, theMidgetRGCmosaic.inputConeMosaic, []);

        % Save data
        fName = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
            ZernikeDataBase, newSubjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2));

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
            mosaicSizeDegs = 1.2;
        elseif (radialEccDegs <= 10)
            mosaicSizeDegs = 1.4;
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
