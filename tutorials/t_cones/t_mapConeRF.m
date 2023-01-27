% Demo usage of RF mapping using the subspace RF mapping method 
%
% Description:
%    Shows how to use the rfMappingStimulusGenerator package to map 
%    the visual spatial RFs of cones using the subspace RF mapping method.
%

% History:
%    09/28/22  NPC  ISETBIO Team, Copyright 2022 Wrote it.

function t_mapConeRF
    
    mosaicEccDegs = [0 -0];
    mosaicSizeDegs = 0.3*[1 1];
    stimSizeDegs = 0.3;

    whichEye = 'right eye';

    % Set cone aperture modifiers
    sigmaGaussian = 0.204;
    coneApertureModifiers = struct(...
            'smoothLocalVariations', true, ...
            'sigma',  sigmaGaussian, ...
            'shape', 'Gaussian');

    % Generate a @cMosaic object
    theConeMosaic = cMosaic(...
        'sourceLatticeSizeDegs', 60, ...
        'eccentricityDegs', mosaicEccDegs, ...
        'sizeDegs', mosaicSizeDegs, ...
        'whichEye', whichEye, ...
        'overlappingConeFractionForElimination', 0.01, ...
        'rodIntrusionAdjustedConeAperture', true, ...
        'coneApertureModifiers', coneApertureModifiers ...
        );

    % Optical params
    mosaicCenterPositionDegs = mean(theConeMosaic.coneRFpositionsDegs,1);
    opticsParams = struct(...
        'positionDegs',mosaicCenterPositionDegs, ...  % (x,y) eccentricity for the PSF, in degrees
        'ZernikeDataBase', 'Polans2015', ...
        'examinedSubjectRankOrder', 10, ...
        'refractiveErrorDiopters', 0.0, ...   % use -999 for optics that do not subtract the central refraction
        'analyzedEye', whichEye, ...
        'subjectRankingEye', 'right eye', ...
        'pupilDiameterMM', 3.0, ...
        'wavefrontSpatialSamples', 501, ...
        'psfUpsampleFactor', 1 ...
        );

    % Generate optics given the optical params
    [thePSFData, ~, ~, ~, theOI] = RetinaToVisualFieldTransformer.computeVlambdaWeightedPSF(...
            opticsParams, theConeMosaic, theConeMosaic.wave);

    % Generate a presentation display with a desired resolution
    retinalImageResolutionDegs = thePSFData.psfSupportXdegs(2)-thePSFData.psfSupportXdegs(1);
    viewingDistanceMeters = 4;
    theDisplay = rfMappingStimulusGenerator.presentationDisplay(...
        theConeMosaic.wave, retinalImageResolutionDegs, viewingDistanceMeters);

    % Stim params for the RF mapping
    stimParams = struct(...
        'backgroundLuminanceCdM2', 50.0, ...
        'backgroundChromaticity', [0.301 0.301], ...
        'coneContrasts', [1 1 1], ...
        'contrast', 0.75, ...
        'pixelSizeDegs', retinalImageResolutionDegs, ...
        'stimSizeDegs', stimSizeDegs, ...
        'wavelengthSupport', displayGet(theDisplay, 'wave'), ...
        'viewingDistanceMeters', displayGet(theDisplay, 'viewing distance') ...
        );


    % Hartley (RF mapping) spatial patterns
    omega = 9;
    % Compute spatial modulation patterns for the Hartley set
    HartleySpatialModulationPatterns = ...
        rfMappingStimulusGenerator.HartleyModulationPatterns(...
        omega, stimParams.stimSizeDegs, stimParams.pixelSizeDegs);
    

    % Generate scenes for the Hartley patterns
    [theRFMappingStimulusScenes, theNullStimulusScene, spatialSupportDegs] = rfMappingStimulusGenerator.generateStimulusFramesOnPresentationDisplay(...
        theDisplay, stimParams, HartleySpatialModulationPatterns, ...
        'validateScenes', false);

    % Generate scenes for the inverse-polarity Hartley patterns
    theReversePolarityRFMappingStimulusScenes = rfMappingStimulusGenerator.generateStimulusFramesOnPresentationDisplay(...
        theDisplay, stimParams, -HartleySpatialModulationPatterns, ...
        'validateScenes', false);

    % Preallocate memory
    conesNum = numel(theConeMosaic.coneTypes);
    nStim = numel(theRFMappingStimulusScenes);
    pixelsNum = size(HartleySpatialModulationPatterns,2);
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

    % Visualize the cone mosaic, the employed stimuli, the PSF and the
    % visually mapped cone RF maps
    visualizeMosaicStimuliAndMappedRFs(theConeMosaic, theConeMosaicExcitation, thePSFData, theDisplay, spatialSupportDegs, theRFMappingStimulusScenes, theRFmaps, maxRFconeType, stimParams);

end

function visualizeMosaicStimuliAndMappedRFs(theConeMosaic, theConeMosaicExcitation, thePSFData, theDisplay, spatialSupportDegs, theRFMappingStimulusScenes, theRFmaps, maxRFconeType, stimParams)
    mosaicCenterPositionDegs = mean(theConeMosaic.coneRFpositionsDegs,1);
    xLims = mosaicCenterPositionDegs(1) + round(0.4*stimParams.stimSizeDegs*[-1 1]*100)/100;
    yLims = mosaicCenterPositionDegs(2) + round(0.4*stimParams.stimSizeDegs*[-1 1]*100)/100;

    xTicks = mosaicCenterPositionDegs(1) + round(0.4*stimParams.stimSizeDegs*[-1:0.5:1]*100)/100;
    yTicks = mosaicCenterPositionDegs(2) + round(0.4*stimParams.stimSizeDegs*[-1:0.5:1]*100)/100;


    videoOBJ = VideoWriter('SubSpaceRFmapping', 'MPEG-4');
    videoOBJ.FrameRate = 10;
    videoOBJ.Quality = 100;
    videoOBJ.open();

    hFig = figure(1);clf;
    set(hFig, 'Resize', 'off', 'Position', [300 500 1270 1000], 'Color', [1 1 1]);
    fontSize = 14;

    axMosaic = subplot(2,3,1);
    theConeMosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', axMosaic, ...
        'domain', 'degrees', ...
        'domainVisualizationLimits', [xLims(1) xLims(2) yLims(1) yLims(2)], ...
        'domainVisualizationTicks', struct('x', xTicks, 'y', yTicks), ...
        'visualizedConeAperture', 'lightCollectingArea4sigma', ...
        'fontSize', fontSize ...
        );


    ax = subplot(2,3,2);
    imagesc(ax, thePSFData.psfSupportXdegs, thePSFData.psfSupportYdegs, thePSFData.vLambdaWeightedData);
    hold(ax, 'on');
    plot(ax, [thePSFData.psfSupportXdegs(1) thePSFData.psfSupportXdegs(end)], [0 0], 'g-');
    plot(ax, [0 0], [thePSFData.psfSupportYdegs(1) thePSFData.psfSupportYdegs(end)], 'g-');
    hold(ax, 'off');
    axis(ax, 'image'); axis(ax, 'xy');
    set(ax, 'XLim', xLims-mosaicCenterPositionDegs(1), 'YLim', yLims-mosaicCenterPositionDegs(2), ...
        'XTick', xTicks-mosaicCenterPositionDegs(1), 'YTick', yTicks-mosaicCenterPositionDegs(2), ...
        'Color', [0 0 0]);
    set(ax, 'FontSize', fontSize);
    colormap(ax, brewermap(1024, 'greys'));
    xlabel(ax, 'space (degrees)');
    title(ax, 'PSF (vLambda-weighted)');

    % Plot the stimuli
    halfLutEntriesNum = 500; 
    cMap = brewermap(2*halfLutEntriesNum+1, '*greys');
    ax1 = subplot(2,3,5);
    ax2 = subplot(2,3,4);


    for iFrame = 1:numel(theRFMappingStimulusScenes)
       % [~, sceneSRGBimage] = ...
        %        sceneRepresentations(theRFMappingStimulusScenes{iFrame}, theDisplay);
        sceneSRGBimage = sceneGet(theRFMappingStimulusScenes{iFrame}, 'RGBimage');
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
        
        % Anatomical cone aperture
        R = sqrt(XX.^2 + YY.^2);

        % Gaussian cone aperture
        theConeAperture = exp(-(R/Rc).^2);

        % Pillbox cone aperture
        %theConeAperture = R * 0;
        %theConeAperture(R <= 0.5*theConeMosaic.coneApertureDiametersDegs(iCone)) = 1;

        % The visual cone aperture
        theVisualConeAperture = conv2(theConeAperture, thePSFData.vLambdaWeightedData, 'same');

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
        legend(ax5, {'anatomical cone aperture', 'conv(cone aperture, vLambda PSF)', 'visual cone RF map'}, 'Location', 'North');
        axis(ax5, 'square');
        grid(ax5, 'on'); box(ax5, 'on');
        set(ax5, 'XLim', xLims, 'YLim', [-0.02 1.4], 'YTick', 0:0.1:1.4, 'XTick', xTicks, 'XTickLabels', sprintf('%2.2f\n', xTicks), 'FontSize', 15);
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
        hold(ax4, 'on');
        plot(ax4, xCross*[1 1], yCross + Rc *[-1 1], 'b-', 'LineWidth', 1.0);
        plot(ax4, xCross + Rc*[-1 1], yCross*[1 1], 'b-', 'LineWidth', 1.0);
        hold(ax4, 'off');
        set(ax4, 'CLim', [-1 1], 'XLim', xLims, 'YLim', yLims, ...
            'Color', squeeze(cMap(halfLutEntriesNum+1,:)));
        axis(ax4, 'image'); axis(ax4, 'xy');
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


function [sceneLRGBimage, sceneSRGBimage, sceneLMScontrastsImage, sceneLMSexcitationsImage] = ...
    sceneRepresentations(theScene, presentationDisplay)

    emittedRadianceImage = sceneGet(theScene, 'energy');
    displaySPDs = displayGet(presentationDisplay, 'spd');
    displayWavelengths = displayGet(presentationDisplay, 'wave');
    [sceneLRGBimage, sceneSRGBimage] = displayRadianceToDisplayRGB(emittedRadianceImage, displaySPDs);
        
    % Load the 2-deg Stockman cone fundamentals on a wavelength support matching the display
    coneFundamentals = ieReadSpectra(fullfile(isetbioDataPath,'human','stockman'), displayWavelengths);

    % Compute the LMS cone contrasts of the emitted radiance image
    [sceneLMScontrastsImage, sceneLMSexcitationsImage] = displayRadianceToLMS(emittedRadianceImage, coneFundamentals);
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