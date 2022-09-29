function t_mapConeRF

    
    mosaicEccDegs = [0.5 -1.2];
    mosaicSizeDegs = [0.1 0.1];
    stimSizeDegs = 0.15;

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
        'wavefrontSpatialSamples', 501, ...
        'psfUpsampleFactor', 1 ...
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
    omega = 11;
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
            r = theConeMosaicExcitation(iFrame,iCone) - theReversePolarityConeMosaicExcitation(iCone);
            theRFmaps(iCone,:,:) = theRFmaps(iCone,:,:) + ...
                HartleySpatialModulationPatterns(iFrame,:,:) * r;
        end
    end
    theRFmaps = 1/(2*nStim)*theRFmaps;

    % Compute maxRF for each cone type
    maxRFconeType = zeros(1,4);
    for iCone = 1:conesNum
         maxRFconeType(theConeMosaic.coneTypes(iCone)) = max([...
             maxRFconeType(theConeMosaic.coneTypes(iCone))  ...
             max(max(squeeze(abs(theRFmaps(iCone,:,:)))))]);
    end

    visualizeMosaicPSFandStimuli(theConeMosaic, theConeMosaicExcitation, HartleySpatialModulationPatterns, thePSFData, theDisplay, spatialSupportDegs, theRFMappingStimulusScenes, theRFmaps, maxRFconeType, stimParams);

end

function visualizeMosaicPSFandStimuli(theConeMosaic, theConeMosaicExcitation, HartleySpatialModulationPatterns, thePSFData, theDisplay, spatialSupportDegs, theRFMappingStimulusScenes, theRFmaps, maxRFconeType, stimParams)
    mosaicCenterPositionDegs = mean(theConeMosaic.coneRFpositionsDegs,1);
    xLims = mosaicCenterPositionDegs(1) + round(0.4*stimParams.stimSizeDegs*[-1 1]*100)/100;
    yLims = mosaicCenterPositionDegs(2) + round(0.4*stimParams.stimSizeDegs*[-1 1]*100)/100;

    xTicks = mosaicCenterPositionDegs(1) + round(0.4*stimParams.stimSizeDegs*[-1:0.5:1]*100)/100;
    yTicks = mosaicCenterPositionDegs(2) + round(0.4*stimParams.stimSizeDegs*[-1:0.5:1]*100)/100;


    videoOBJ = VideoWriter('SubSpaceRFmapping', 'MPEG-4');
    videoOBJ.FrameRate = 10;
    videoOBJ.Quality = 100;
    videoOBJ.open();

    hFig = figure(1); clf;
    fontSize = 15;

    axMosaic = subplot(2,3,1);
    theConeMosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', axMosaic, ...
        'domain', 'degrees', ...
        'domainVisualizationLimits', [xLims(1) xLims(2) yLims(1) yLims(2)], ...
        'domainVisualizationTicks', struct('x', xTicks, 'y', yTicks), ...
        'visualizedConeAperture', 'lightCollectingArea4sigma', ...
        'fontSize', fontSize, ...
        'noXLabel', true ...
        );


    ax = subplot(2,3,2);
    imagesc(ax, thePSFData.supportXdegs, thePSFData.supportYdegs, thePSFData.data);
    axis(ax, 'image');
    set(ax, 'XLim', xLims-mosaicCenterPositionDegs(1), 'YLim', yLims-mosaicCenterPositionDegs(2), ...
        'XTick', xTicks-mosaicCenterPositionDegs(1), 'YTick', yTicks-mosaicCenterPositionDegs(2), ...
        'Color', [0 0 0]);
    set(ax, 'FontSize', fontSize);
    colormap(ax, 'gray');
    xlabel(ax, 'space (degrees)');
    title(ax, 'PSF (vLambda-weighted)');

    % Plot the stimuli
    halfLutEntriesNum = 500; cMap = brewermap(2*halfLutEntriesNum+1, '*greys');
    ax1 = subplot(2,3,5);
    ax2 = subplot(2,3,4);


    for iFrame = 1:numel(theRFMappingStimulusScenes)
        [~, sceneSRGBimage] = ...
                sceneRepresentations(theRFMappingStimulusScenes{iFrame}, theDisplay);

        image(ax1, spatialSupportDegs+mosaicCenterPositionDegs(1), spatialSupportDegs+mosaicCenterPositionDegs(2), sceneSRGBimage);
        axis(ax1, 'image');
        set(ax1, 'XLim', xLims,  'YLim', yLims, 'XTick', [], 'YTick', [], 'Color', sceneSRGBimage(1,1,:));
        set(ax1, 'FontSize', fontSize);
        title(ax1, sprintf('stimulus scene %d/%d',iFrame, numel(theRFMappingStimulusScenes)));

        theConeMosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax2, ...
            'domain', 'degrees', ...
            'domainVisualizationLimits', [xLims(1) xLims(2) yLims(1) yLims(2)], ...
            'domainVisualizationTicks', struct('x', xTicks, 'y', yTicks), ...
            'visualizedConeAperture', 'lightCollectingArea4sigma', ...
            'activation', squeeze(theConeMosaicExcitation(iFrame,:)), ...
            'plotTitle', sprintf('response to stim %d/%d (4 sigma)',iFrame, numel(theRFMappingStimulusScenes)), ...
            'fontSize', fontSize ...
        );

        drawnow;
        videoOBJ.writeVideo(getframe(hFig));
           
    end

    

    % Plot the RFs
    maxRF = max(abs(theRFmaps(:)));
    conesNum = size(theRFmaps,1);
    halfLutEntriesNum = 500; cMap = brewermap(2*halfLutEntriesNum+1, '*RdBu');


    [X,Y] = meshgrid(spatialSupportDegs,spatialSupportDegs);
    R = sqrt(X.^2 + Y.^2);

    ax4 = subplot(2,3,6);
    ax5 = subplot(2,3,3);

    for iCone = 1:conesNum

        xCross = theConeMosaic.coneRFpositionsDegs(iCone,1);
        yCross = theConeMosaic.coneRFpositionsDegs(iCone,2);
        Rc = theConeMosaic.coneApertureDiametersDegs(iCone) * theConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor;

        [XX,YY] = meshgrid(spatialSupportDegs+mosaicCenterPositionDegs(1)-xCross,spatialSupportDegs+mosaicCenterPositionDegs(2)-yCross);
        R = sqrt(XX.^2 + YY.^2);
        %theConeAperture = exp(-(R/Rc).^2);
        theConeAperture = R * 0;
        theConeAperture(R <= 0.5*theConeMosaic.coneApertureDiametersDegs(iCone)) = 1;

        theVisualConeAperture = conv2(theConeAperture, thePSFData.data, 'same');

        [~,maxIndex] = max(theVisualConeAperture(:));
        [row,col] = ind2sub(size(R), maxIndex);
        cla(ax5);
        shadedAreaPlot(ax5,spatialSupportDegs+mosaicCenterPositionDegs(1), theConeAperture(row,:)/max(theConeAperture(:)), 0, ...
            [1 1 0.5], [0 0 0], 0.7, 1.0, '-');
        hold(ax5, 'on');
        plot(ax5, spatialSupportDegs+mosaicCenterPositionDegs(1), theVisualConeAperture(row,:)/max(theVisualConeAperture(:)), 'r-', 'LineWidth', 1.5);
        
        theCurrentRFmap = squeeze(theRFmaps(iCone,:,:));
        plot(ax5, spatialSupportDegs+mosaicCenterPositionDegs(1), theCurrentRFmap(row,:)/max(theCurrentRFmap(:)), 'b-', 'LineWidth', 1.5);
        hold(ax5, 'off');
        legend(ax5, {'cone aperture', 'conv(cone aperture, vLambda PSF)', 'measured cone RF map'}, 'Location', 'NorthOutside');
        axis(ax5, 'square');
        grid(ax5, 'on');
        set(ax5, 'XLim', xLims, 'YLim', [-0.1 1.1], 'YTick', 0:0.1:1, 'XTick', xTicks, 'XTickLabels', sprintf('%2.2f\n', xTicks), 'FontSize', 15);
        xlabel(ax5, 'space (degrees)');

        theConeMosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', axMosaic, ...
            'domain', 'degrees', ...
            'domainVisualizationLimits', [xLims(1) xLims(2) yLims(1) yLims(2)], ...
            'domainVisualizationTicks', struct('x', xTicks, 'y', yTicks), ...
            'visualizedConeAperture', 'lightCollectingArea4sigma', ...
            'labelConesWithIndices', iCone, ...
            'fontSize', fontSize ...
        );
        

        zLevels = [-0.95:0.1:-0.05 0.05:0.1:0.95];
        contourf(ax4, spatialSupportDegs+mosaicCenterPositionDegs(1), ...
                     spatialSupportDegs+mosaicCenterPositionDegs(2), ...
                     squeeze(theRFmaps(iCone,:,:))/maxRFconeType(theConeMosaic.coneTypes(iCone)), zLevels);
        axis(ax4, 'image');
        hold(ax4, 'on');
        
        
        plot(ax4, xCross*[1 1], yCross + Rc *[-1 1], 'b-', 'LineWidth', 1.0);
        plot(ax4, xCross + Rc*[-1 1], yCross*[1 1], 'b-', 'LineWidth', 1.0);
        hold(ax4, 'off');
        set(ax4, 'CLim', [-1 1], 'XLim', xLims, 'YLim', yLims, ...
            'Color', squeeze(cMap(halfLutEntriesNum+1,:)));
        axis(ax4, 'xy');
        set(ax4, 'XTick', xTicks, 'YTick', yTicks, 'XTickLabels', sprintf('%2.2f\n', xTicks), 'YTickLabels', sprintf('%2.2f\n', yTicks));
        grid(ax4, 'on');
        set(ax4, 'FontSize', fontSize);
        xlabel(ax4, 'space (degrees)');
        title(ax4, sprintf('visual RF map - cone %d/%d',iCone, conesNum));
        colormap(ax4, cMap);


        drawnow;
        videoOBJ.writeVideo(getframe(hFig));
        
    end
    videoOBJ.close();
end

function shadedAreaPlot(ax,x,y, baseline, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)
    x = [x fliplr(x)];
    y = [y y*0+baseline];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, ...
        'FaceAlpha', faceAlpha, 'LineWidth', lineWidth, 'LineStyle', lineStyle);
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
            sIndex = sIndex + 1;
            H(sIndex,:,:) = f;
            lIndices(sIndex) = lIndex-omega;
            mIndices(sIndex) = mIndex-omega;
        end
    end
    H = H / max(H(:));
end





