function t_mapConeRF

    
    mosaicEccDegs = [0.5 1.2];
    mosaicSizeDegs = [0.1 0.1];
    stimSizeDegs = 0.1;
    whichEye = 'right eye';

    % Set cone aperture modifiers
    sigmaGaussian = 0.204;
    coneApertureModifiers = struct(...
            'smoothLocalVariations', true, ...
            'sigma',  sigmaGaussian, ...
            'shape', 'Gaussian');

    theConeMosaic = cMosaic(...
        'sourceLatticeSizeDegs', 60, ...
        'eccentricityDegs', mosaicEccDegs, ...
        'sizeDegs', mosaicSizeDegs, ...
        'whichEye', whichEye, ...
        'overlappingConeFractionForElimination', 0.01, ...
        'rodIntrusionAdjustedConeAperture', true, ...
        'coneApertureModifiers', coneApertureModifiers ...
        );

    mosaicCenterPositionDegs = mean(theConeMosaic.coneRFpositionsDegs,1);

    opticsParams = struct(...
        'positionDegs',mosaicCenterPositionDegs, ...  % (x,y) eccentricity for the PSF, in degrees
        'ZernikeDataBase', 'Artal2012', ...
        'examinedSubjectRankOrder', 2, ...
        'refractiveErrorDiopters', 0.0, ...   % use -999 for optics that do not subtract the central refraction
        'analyzedEye', whichEye, ...
        'subjectRankingEye', 'right eye', ...
        'pupilDiameterMM', 3.0, ...
        'wavefrontSpatialSamples', 401, ...
        'psfUpsampleFactor', 5 ...
        );

    [thePSFData, ~, ~, theOI] = RetinaToVisualFieldTransformer.computeVlambdaWeightedPSF(...
            opticsParams, theConeMosaic, theConeMosaic.wave);

    retinalImageResolutionDegs = thePSFData.supportXdegs(2)-thePSFData.supportXdegs(1);
    viewingDistanceMeters = 4;
    theDisplay = rfMappingStimulusGenerator.presentationDisplay(theConeMosaic.wave, retinalImageResolutionDegs, viewingDistanceMeters);

    stimParams = struct(...
        'backgroundLuminanceCdM2', 50.0, ...
        'backgroundChromaticity', [0.301 0.301], ...
        'coneContrasts', [1 1 1], ...
        'contrast', 0.75, ...
        'pixelSizeDegs', retinalImageResolutionDegs, ...
        'stimSizeDegs', stimSizeDegs, ...
        'spatialFrequencyCPD', 120, ...
        'orientationDegs', 20, ...
        'spatialPhaseDegs', 0, ...
        'wavelengthSupport', displayGet(theDisplay, 'wave'), ...
        'viewingDistanceMeters', displayGet(theDisplay, 'viewing distance') ...
        );


    % Generate Hartley modulation patterns
    omega = 4;
    % Compute spatial modulation patterns for the Hartley set
    pixelsNum  = round(stimParams.stimSizeDegs / stimParams.pixelSizeDegs);
    pixelIndex = 0:(pixelsNum-1);
    [X,Y] = meshgrid(pixelIndex/pixelsNum, pixelIndex/pixelsNum);
    [HartleySpatialModulationPatterns, lIndices, mIndices] = HartleyModulationPatterns(X,Y,omega);

    % Generate scenes for the Hartley patterns
    [theRFMappingStimulusScenes, theNullStimulusScene, spatialSupportDegs] = rfMappingStimulusGenerator.generateStimulusMappingFramesOnPresentationDisplay(...
        theDisplay, stimParams, HartleySpatialModulationPatterns, ...
        'validateScenes', false);
    % Generate the scenes for the reverse Hartley patterns
    theReversePolarityRFMappingStimulusScenes = rfMappingStimulusGenerator.generateStimulusMappingFramesOnPresentationDisplay(...
        theDisplay, stimParams, -HartleySpatialModulationPatterns, ...
        'validateScenes', false);

    % Preallocate memory
    conesNum = numel(theConeMosaic.coneTypes);
    nStim = numel(theRFMappingStimulusScenes);
    theRFmaps = zeros(conesNum, pixelsNum, pixelsNum);
    theConeMosaicExcitation = zeros(nStim, conesNum);

    % Compute the cone mosaic responses and build - up the RF
    for iFrame = 0:nStim
        if (iFrame == 0)
            theOI = oiCompute(theNullStimulusScene, theOI);
            theNullStimConeMosaicExcitation = squeeze(theConeMosaic.compute(theOI));
            continue;
        end

        fprintf('Computing mosaic response to stim %d of %d\n', iFrame, nStim);
        % The forward polarity responses
        theOI = oiCompute(theRFMappingStimulusScenes{iFrame}, theOI);
        theConeMosaicExcitation(iFrame,:) = (squeeze(theConeMosaic.compute(theOI)) - theNullStimConeMosaicExcitation)./theNullStimConeMosaicExcitation;
        
        % The reverse polarity responses
        theOI = oiCompute(theReversePolarityRFMappingStimulusScenes{iFrame}, theOI);
        theReversePolarityConeMosaicExcitation = (squeeze(theConeMosaic.compute(theOI)) - theNullStimConeMosaicExcitation)./theNullStimConeMosaicExcitation;

        % Update the RF map for each cone
        parfor iCone = 1:conesNum
            theRFmaps(iCone,:,:) = theRFmaps(iCone,:,:) + ...
                HartleySpatialModulationPatterns(iFrame,:,:) * (theConeMosaicExcitation(iFrame,iCone) - theReversePolarityConeMosaicExcitation(iCone));
        end
    end

    % Compute RF for each cone
    theRFmaps = 1/(2*nStim)*theRFmaps;

    for iCone = 1:conesNum
        maxRF(iCone) = max(max(squeeze(abs(theRFmaps(iCone,:,:)))));
        fprintf('max RF map (cone %d): %f\n', iCone, maxRF(iCone));
    end

    visualizeMosaicPSFandStimuli(theConeMosaic, theConeMosaicExcitation, HartleySpatialModulationPatterns, thePSFData, theDisplay, spatialSupportDegs, theRFMappingStimulusScenes, theRFmaps);

end

function visualizeMosaicPSFandStimuli(theConeMosaic, theConeMosaicExcitation, HartleySpatialModulationPatterns, thePSFData, theDisplay, spatialSupportDegs, theRFMappingStimulusScenes, theRFmaps)
    mosaicCenterPositionDegs = mean(theConeMosaic.coneRFpositionsDegs,1);
    xLims = mosaicCenterPositionDegs(1) + theConeMosaic.sizeDegs(1)*0.5*[-1 1];
    yLims = mosaicCenterPositionDegs(2) + theConeMosaic.sizeDegs(2)*0.5*[-1 1];

    hFig = figure(1); clf;
    fontSize = 15;

    ax0 = subplot(2,4,1);
    theConeMosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'domain', 'degrees', ...
        'domainVisualizationLimits', [xLims(1) xLims(2) yLims(1) yLims(2)], ...
        'visualizedConeAperture', 'geometricArea', ...
        'plotTitle', 'geometricArea', ...
        'fontSize', fontSize ...
        );

    ax = subplot(2,4,2);
    theConeMosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'domain', 'degrees', ...
        'domainVisualizationLimits', [xLims(1) xLims(2) yLims(1) yLims(2)], ...
        'visualizedConeAperture', 'coneSpacing', ...
        'plotTitle', 'coneSpacing', ...
        'fontSize', fontSize ...
        );

    ax = subplot(2,4,3);
    imagesc(ax, thePSFData.supportXdegs, thePSFData.supportYdegs, thePSFData.data);
    axis(ax, 'image');
    set(ax, 'XLim', xLims-mosaicCenterPositionDegs(1), 'YLim', yLims-mosaicCenterPositionDegs(2), 'Color', [0 0 0]);
    set(ax, 'FontSize', fontSize);
    title(ax, 'PSF');

    % Plot the stimuli
    halfLutEntriesNum = 500; cMap = brewermap(2*halfLutEntriesNum+1, '*greys');
    ax1 = subplot(2,4,5);
    ax2 = subplot(2,4,6);
    ax3 = subplot(2,4,7);
    for iFrame = 1:numel(theRFMappingStimulusScenes)
        [~, sceneSRGBimage] = ...
                sceneRepresentations(theRFMappingStimulusScenes{iFrame}, theDisplay);

        image(ax1, spatialSupportDegs+mosaicCenterPositionDegs(1), spatialSupportDegs+mosaicCenterPositionDegs(2), sceneSRGBimage);
        axis(ax1, 'image'); axis(ax1, 'ij');
        set(ax1, 'XLim', xLims,  'YLim', yLims, 'Color', sceneSRGBimage(1,1,:));
        set(ax1, 'FontSize', fontSize);
        title(ax1, sprintf('stimulus scene %d/%d',iFrame, numel(theRFMappingStimulusScenes)));

        theConeMosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax2, ...
            'domain', 'degrees', ...
            'domainVisualizationLimits', [xLims(1) xLims(2) yLims(1) yLims(2)], ...
            'visualizedConeAperture', 'lightCollectingArea4sigma', ...
            'activation', squeeze(theConeMosaicExcitation(iFrame,:)), ...
            'plotTitle', sprintf('response to stim %d/%d (4 sigma)',iFrame, numel(theRFMappingStimulusScenes)), ...
            'fontSize', fontSize ...
        );

        imagesc(ax3, spatialSupportDegs+mosaicCenterPositionDegs(1), spatialSupportDegs+mosaicCenterPositionDegs(2), ...
            squeeze(HartleySpatialModulationPatterns(iFrame,:,:)));
        axis(ax3, 'image'); axis(ax3, 'ij');
        set(ax3, 'XLim', xLims, 'YLim', yLims, 'CLim', [-1 1], 'Color', squeeze(cMap(halfLutEntriesNum+1,:)));
        set(ax3, 'FontSize', fontSize);
        title(ax3, sprintf('stimulus spatial modulation %d/%d',iFrame, numel(theRFMappingStimulusScenes)));
        colormap(ax3, cMap);
        drawnow;
        pause
    end

    

    % Plot the RFs
    maxRF = max(abs(theRFmaps(:)));
    conesNum = size(theRFmaps,1);
    halfLutEntriesNum = 500; cMap = brewermap(2*halfLutEntriesNum+1, '*RdBu');

    ax4 = subplot(2,4,8);
    for iCone = 1:conesNum

        theConeMosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax0, ...
            'domain', 'degrees', ...
            'domainVisualizationLimits', [xLims(1) xLims(2) yLims(1) yLims(2)], ...
            'visualizedConeAperture', 'lightCollectingArea4sigma', ...
            'activation', squeeze(theConeMosaicExcitation(iFrame,:)), ...
            'labelConesWithIndices', iCone, ...
            'plotTitle', sprintf('response to stim %d/%d (4 sigma)',iFrame, numel(theRFMappingStimulusScenes)), ...
            'fontSize', fontSize ...
        );
        

        imagesc(ax4, spatialSupportDegs+mosaicCenterPositionDegs(1), spatialSupportDegs+mosaicCenterPositionDegs(2), squeeze(theRFmaps(iCone,:,:))/maxRF);
        hold(ax, 'on');
        
        xCross = theConeMosaic.coneRFpositionsDegs(iCone,1);
        yCross = theConeMosaic.coneRFpositionsDegs(iCone,2);
        Rc = theConeMosaic.coneApertureDiametersDegs(iCone) * theConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor;
        plot(ax4, xCross*[1 1], yCross + Rc *[-1 1], 'k-', 'LineWidth', 1.0);
        plot(ax4, xCross + Rc*[-1 1], yCross*[1 1], 'k-', 'LineWidth', 1.0);
        hold(ax4, 'off');
        axis(ax4, 'image'); axis(ax, 'ij');

        set(ax4, 'CLim', [-1 1], 'XLim', xLims, 'YLim', yLims, 'Color', squeeze(cMap(halfLutEntriesNum+1,:)));
        set(ax4, 'FontSize', fontSize);
        title(ax4, sprintf('cone %d/%d',iCone, conesNum));
        colormap(ax4, cMap);


        drawnow;
        pause
    end

end


function [H, lIndices, mIndices] = HartleyModulationPatterns(X,Y,omega)

    H = zeros((2*omega+1)^2, size(X,1), size(X,2));
    lIndices = zeros((2*omega+1)^2,1);
    mIndices = lIndices;
    sIndex = 0;
    for mIndex = 0:(2*omega)
        for lIndex = 0:(2*omega)
            a = 2*pi*((lIndex-omega)*X + (mIndex-omega)*Y);
            f = sin(a)+cos(a);
            f(1:80,:) = 0;
            sIndex = sIndex + 1;
            H(sIndex,:,:) = f;
            lIndices(sIndex) = lIndex-omega;
            mIndices(sIndex) = mIndex-omega;
        end
    end
    H = H / max(H(:));
end





