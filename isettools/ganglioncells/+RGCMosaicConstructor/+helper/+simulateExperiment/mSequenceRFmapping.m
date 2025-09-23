function [theMinorSigmasDegs, theMajorSigmasDegs, theTemporalEquivalentEccentricitiesDegs, ...
          theMinorSigmasMicrons, theMajorSigmasMicrons, theTemporalEquivalentEccentricitiesMMs] = mSequenceRFmapping(...
    mRGCMosaicParamsStruct, opticsParamsStruct, RFmappingParamsStruct, ...
    performComputeInputConeMosaicResponsesAction, ...
    performComputeMRGCMosaicResponsesAction, ...
    performVisualizeMRGCMosaicRFmaps, ...
    varargin)
    
    % Parse input
    p = inputParser;
    % Optional params
    p.addParameter('visualizeInputConeMosaicResponses', false, @islogical);
    p.parse(varargin{:});
    visualizeInputConeMosaicResponses = p.Results.visualizeInputConeMosaicResponses;

    theMinorSigmas = [];
    theMajorSigmas = [];

    % Generate filepath for center-connected mosaic
    if (mRGCMosaicParamsStruct.rfCentersOverlap)
        centerConnectivityStageDescriptor = 'center connected with overlap';
        surroundConnectivityStageDescriptor = 'surround connected with center overlap';
    else
        centerConnectivityStageDescriptor = 'center connected';
        surroundConnectivityStageDescriptor = 'surround connected';
    end
    % Generate mosaic filename
    [theMRGCMosaicFullFileName, ~, theMRGCMosaicFileName, mRGCMosaicSubDir] = ...
        RGCMosaicConstructor.filepathFor.exportedMosaicFileName(mRGCMosaicParamsStruct, centerConnectivityStageDescriptor);

    % Assemble MsequenceRFmappingResponses  filename: add sub-directory- inputConeMosaicMsequenceRFmappingResponses
    theInputConeMosaicMsequenceRFmappingResponsesFullFileName = strrep(theMRGCMosaicFullFileName, 'SLIM', 'SLIM/inputConeMosaicMsequenceRFmappingResponses');

    % Update MsequenceRFmappingResponses filename to include the optics employed & stimulus position
    opticsString = sprintf('OpticsAt_%2.1f_%2.1f_%s-%d_%s', ...
                RFmappingParamsStruct.positionDegs(1), ...
                RFmappingParamsStruct.positionDegs(2), ...
                opticsParamsStruct.ZernikeDataBase, ...
                opticsParamsStruct.subjectID, ...
                opticsParamsStruct.modification);
    theInputConeMosaicMsequenceRFmappingResponsesFullFileName = strrep(theInputConeMosaicMsequenceRFmappingResponsesFullFileName, ...
        '.mat', sprintf('_%s.mat', opticsString));

    % Update MsequenceRFmappingResponses filename to include the stimulus chromaticity info 
    chromaParamsStruct = struct(...
        'stimulusChromaticity', RFmappingParamsStruct.chromaticity, ...
        'coneFundamentalsOptimizedForStimPosition', RFmappingParamsStruct.coneFundamentalsOptimizedForStimPosition);
    [~, ~, theInputConeMosaicMsequenceRFmappingResponsesFullFileName] = ...
        RGCMosaicConstructor.helper.queryUserFor.stimulusChromaticity(theInputConeMosaicMsequenceRFmappingResponsesFullFileName, ...
        'doNotQueryUserInsteadUseThisChromaParams', chromaParamsStruct);

    % Generate corresponding mRGCMosaic STF responses filaname
    theMRGCMosaicMsequenceRFmappingResponsesFullFileName = ...
        strrep(theInputConeMosaicMsequenceRFmappingResponsesFullFileName, 'inputConeMosaicMsequenceRFmappingResponses', 'MsequenceRFmappingResponses');

    % Remove the connectiviStageDescriptor from  theInputConeMosaicMsequenceRFmappingResponsesFullFileName
    theInputConeMosaicMsequenceRFmappingResponsesFullFileName = ...
        strrep(theInputConeMosaicMsequenceRFmappingResponsesFullFileName, mRGCMosaicSubDir, '');

    % Compute input cone mosaic responses
    if (performComputeInputConeMosaicResponsesAction)
        % Load the mosaic so we can crop it (to the size mapped)
        theImportedData = load(theMRGCMosaicFullFileName, 'theMRGCMosaic');

        % Compute extra support for input cone mosaic at the optimization position
        extraSupportDegsForInputConeMosaic = mRGCMosaic.extraSupportDegsForMidgetRGCSurrounds(...
            RFmappingParamsStruct.positionDegs, RFmappingParamsStruct.sizeDegs);

        % Crop the mosaic
        theMRGCMosaic = theImportedData.theMRGCMosaic.cropToSizeAtEccentricity(...
            RFmappingParamsStruct.sizeDegs, RFmappingParamsStruct.positionDegs, ...
            'extraSupportDegsForInputConeMosaic', 0.5*extraSupportDegsForInputConeMosaic, ...
            'visualizeSpatialRelationshipToSourceMosaic', ~true);
        clear 'theImportedData';

        % Generate optics
        [theOI, thePSF] = RGCMosaicConstructor.helper.optics.generate(...
            theMRGCMosaic.inputConeMosaic, ...
            RFmappingParamsStruct.positionDegs, ...
            opticsParamsStruct, ...
            'visualizePSF', true);

        % Determine the stimulus pixel resolution to be a fraction of the minimum cone aperture or cone spacing in the mosaic
        % here, half of the cone spacing
        theMetric = 'cone aperture';  % choose from {'cone aperture' or cone spacing'}
        theFraction = 0.25;

        targetRGCindices =  1:theMRGCMosaic.rgcsNum;
        RFmappingParamsStruct.optimalResolutionDegs = RGCMosaicConstructor.helper.simulateExperiment.stimulusResolutionFromConeApertureOrConeSpacing(...
                theMRGCMosaic, targetRGCindices, theFraction, theMetric);

        % Cover all of the input cone mosaic
        RFmappingParamsStruct.sizeDegs =  max(theMRGCMosaic.inputConeMosaic.sizeDegs);

        % Runc the m-sequence RF mapping
        mSequenceRFmappingCore(...
            theMRGCMosaic, theOI, RFmappingParamsStruct, ...
            theInputConeMosaicMsequenceRFmappingResponsesFullFileName, ...
            theMRGCMosaicMsequenceRFmappingResponsesFullFileName, ...
            'computeInputConeMosaicResponses', true, ...
            'visualizeInputConeMosaicResponses', visualizeInputConeMosaicResponses, ...
            'computeMRGCMosaicResponses', false);
    end % if (performComputeInputConeMosaicResponsesAction)

    if (performComputeMRGCMosaicResponsesAction)
        mSequenceRFmappingCore(...
            [], [], [], ...
            theInputConeMosaicMsequenceRFmappingResponsesFullFileName, ...
            theMRGCMosaicMsequenceRFmappingResponsesFullFileName, ...
            'computeInputConeMosaicResponses', false, ...
            'visualizeInputConeMosaicResponses', false, ...
            'computeMRGCMosaicResponses', true);
    end % if (performComputeInputConeMosaicResponsesAction)

    if (performVisualizeMRGCMosaicRFmaps)
        [theMinorSigmasDegs, theMajorSigmasDegs, theTemporalEquivalentEccentricitiesDegs, ...
        theMinorSigmasMicrons, theMajorSigmasMicrons, theTemporalEquivalentEccentricitiesMMs] = mSequenceRFmappingCore(...
            [], [], [], ...
            theInputConeMosaicMsequenceRFmappingResponsesFullFileName, ...
            theMRGCMosaicMsequenceRFmappingResponsesFullFileName, ...
            'computeInputConeMosaicResponses', false, ...
            'visualizeInputConeMosaicResponses', false, ...
            'computeMRGCMosaicResponses', false, ...
            'visualizeMRGCMosaicRFmaps', true);
    end
end



function [theMinorSigmasDegs, theMajorSigmasDegs, theTemporalEquivalentEccentricitiesDegs, ...
          theMinorSigmasMicrons, theMajorSigmasMicrons, theTemporalEquivalentEccentricitiesMMs] = ...
                mSequenceRFmappingCore(theMRGCMosaic, theOI, ...
                RFmappingParamsStruct, theInputConeMosaicResponsesFullFileName, theMRGCMosaicResponsesFullFileName, varargin)

  p = inputParser;
  p.addParameter('computeInputConeMosaicResponses', false, @islogical);
  p.addParameter('visualizeInputConeMosaicResponses', false, @islogical);
  p.addParameter('computeMRGCMosaicResponses', false, @islogical);
  p.addParameter('visualizeMRGCMosaicRFmaps', false, @islogical);
  % Execute the parser
  p.parse(varargin{:});

  computeInputConeMosaicResponses = p.Results.computeInputConeMosaicResponses;
  computeMRGCMosaicResponses = p.Results.computeMRGCMosaicResponses;
  visualizeInputConeMosaicResponses = p.Results.visualizeInputConeMosaicResponses;
  visualizeMRGCMosaicRFmaps = p.Results.visualizeMRGCMosaicRFmaps;

 if (computeInputConeMosaicResponses)
    inputConeMosaicMsequenceRFmapping(theMRGCMosaic, theOI, ...
        RFmappingParamsStruct, theInputConeMosaicResponsesFullFileName, ...
        visualizeInputConeMosaicResponses);
  end

  if (computeMRGCMosaicResponses)
    mRGCMosaicMsequenceRFmapping(...
        theInputConeMosaicResponsesFullFileName, ...
        theMRGCMosaicResponsesFullFileName, ...
        visualizeInputConeMosaicResponses);
  end

  theMinorSigmasDegs = [];
  theMajorSigmasDegs = [];
  theTemporalEquivalentEccentricitiesDegs = [];
  theMinorSigmasMicrons = [];
  theMajorSigmasMicrons = [];
  theTemporalEquivalentEccentricitiesMMs = [];

  if (visualizeMRGCMosaicRFmaps)
    [theMinorSigmasDegs, theMajorSigmasDegs, theTemporalEquivalentEccentricitiesDegs, ...
    theMinorSigmasMicrons, theMajorSigmasMicrons, theTemporalEquivalentEccentricitiesMMs] = ...
            mRGCMosaicRFmappingVisualizations(theMRGCMosaicResponsesFullFileName);
  end

end

function [theMinorSigmasDegs, theMajorSigmasDegs, theTemporalEquivalentEccentricitiesDegs, ...
         theMinorSigmasMicrons, theMajorSigmasMicrons, theTemporalEquivalentEccentricitiesMMs] = mRGCMosaicRFmappingVisualizations(theMRGCMosaicResponsesFullFileName)

    load(theMRGCMosaicResponsesFullFileName, ...
        'theTemporalEquivalentEccentricitiesDegs', ...
        'theTemporalEquivalentEccentricitiesMMs', ...
        'theRFmaps', ...
        'spatialSupportDegs', ...
        'theFittedGaussianEllipsoids');

    theMinorSigmasDegs = zeros(1,numel(theRFmaps));
    theMajorSigmasDegs = zeros(1,numel(theRFmaps));
    theMinorSigmasMicrons = zeros(1,numel(theRFmaps));
    theMajorSigmasMicrons = zeros(1,numel(theRFmaps));

    hFig = figure(1); clf;
    axIndividual = subplot(1,2,1);
    axEnsemble = subplot(1,2,2);

    XLims = [spatialSupportDegs(1) spatialSupportDegs(end)];
    YLims = [spatialSupportDegs(1) spatialSupportDegs(end)];
    XLims = 0.8 * [-1 1];
    YLims = 0.8 * [-1 1];
    for iRGC = 1:numel(theRFmaps)
        theRawRFmap = theRFmaps{iRGC};
        theFittedGaussianEllipsoid = theFittedGaussianEllipsoids{iRGC};
        theFittedRFmap = theFittedGaussianEllipsoid.rfMap;
        theFittedRFmapSpatialSupportXDegs = theFittedGaussianEllipsoid.xSupportDegs;
        theFittedRFmapSpatialSupportYDegs = theFittedGaussianEllipsoid.ySupportDegs;

        theMinorSigmasDegs(iRGC) = theFittedGaussianEllipsoid.minorSigmaDegs;
        theMajorSigmasDegs(iRGC) = theFittedGaussianEllipsoid.majorSigmaDegs;
        theMinorSigmasMicrons(iRGC) = theFittedGaussianEllipsoid.minorSigmaMicrons;
        theMajorSigmasMicrons(iRGC) = theFittedGaussianEllipsoid.majorSigmaMicrons;

        zLevels = exp(-0.5)*[1 1];
        
        imagesc(axIndividual,spatialSupportDegs, spatialSupportDegs, theRawRFmap/max(theRawRFmap(:)));
        hold(axIndividual, 'on');
        contour(axIndividual, theFittedRFmapSpatialSupportXDegs, theFittedRFmapSpatialSupportYDegs, theFittedRFmap/max(theFittedRFmap(:)), ...
            zLevels, 'Color', [1 0 0], 'LineWidth', 1.5);
        plot(axIndividual, spatialSupportDegs*0, spatialSupportDegs, 'k-');
        plot(axIndividual, spatialSupportDegs, spatialSupportDegs*0, 'k-');
        hold(axIndividual,'off')
        box(axIndividual, 'on');
        
        axis(axIndividual, 'image'); axis(axIndividual, 'xy');
        set(axIndividual, 'CLim', [0 1], 'ZLim', [0 1], 'XLim', XLims, 'YLim', YLims);
        set(axIndividual, 'fontSize', 16);
        xlabel(axIndividual, 'space, x (degs)');
        ylabel(axIndividual, 'space, y (degs)');
        colormap(axIndividual, brewermap(1024, '*greys'));


        contour(axEnsemble, theFittedRFmapSpatialSupportXDegs, theFittedRFmapSpatialSupportYDegs, theFittedRFmap/max(theFittedRFmap(:)), ...
            zLevels, 'Color', [1 0 0], 'LineWidth', 1.0);
        %plot(axEnsemble, spatialSupportDegs*0, spatialSupportDegs, 'k-');
        %plot(axEnsemble, spatialSupportDegs, spatialSupportDegs*0, 'k-');
        hold(axEnsemble,'on')
        axis(axEnsemble, 'image'); axis(axEnsemble, 'xy');
        set(axIndividual, 'XColor', 'none', 'YColor', 'none');
        set(axEnsemble, 'CLim', [0 1], 'ZLim', [0 1], 'XLim', XLims, 'YLim', YLims);
        xlabel(axEnsemble, 'space, x (degs)');
        ylabel(axEnsemble, 'space, y (degs)');
        set(axEnsemble, 'fontSize', 16);
        drawnow;
        pause
    end

    hFig = figure(2); clf;
    sigmaRange = [min(theMinorSigmasDegs) max(theMajorSigmasDegs)];
    plot(theMajorSigmasDegs, theMinorSigmasDegs, 'ro', 'MarkerFaceColor', [1 0.5 0.5]);
    hold on;
    plot(sigmaRange, sigmaRange, 'k-');
    axis 'square'
    set(gca, 'XLim', sigmaRange, 'YLim', sigmaRange);
    xlabel('major axis sigma (degs)')
    ylabel('minor axis sigma (degs)');
end

function mRGCMosaicMsequenceRFmapping(theInputConeMosaicResponsesFullFileName, theMRGCMosaicResponsesFullFileName, visualizeInputConeMosaicResponses)

  fprintf('Loading computed input cone mosaic m-sequence responses from %s.\n', theInputConeMosaicResponsesFullFileName);
  load(theInputConeMosaicResponsesFullFileName, ...
    'theMRGCMosaic', ...
    'RFmappingParamsStruct', ...
    'stimParams', ...
    'spatialSupportDegs', ...
    'mSequenceIndicatorFunctions', ...
    'theInputConeMosaicResponses');
  
    hFig = figure(1000); clf;
    % Prepare figure and axes
    ff = PublicationReadyPlotLib.figureComponents('1x1 giant square mosaic');
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    
    xMin = -15; xMax = -11.0;
    yMin = -2.0; yMax = 2.0;
    domainVisualizationLimits = [xMin xMax yMin yMax];
    domainVisualizationTicks = struct('x', -30:0.2:30, 'y', -20:0.2:20);
    theMRGCMosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', theAxes{1,1}, ...
        'domainVisualizationLimits', domainVisualizationLimits, ...
        'domainVisualizationTicks', domainVisualizationTicks, ...
        'minConeWeightVisualized', exp(-0.5), ...
        'identifiedConeAperture', 'lightCollectingArea4sigma', ... %'lightCollectingAreaCharacteristicDiameter', ...
        'identifiedConeApertureThetaSamples', 30, ...
        'identifyInputCones', ~true, ...
        'identifyPooledCones', ~true, ...
        'inputConesAlpha', 0.9, ...
        'pooledConesLineWidth', 1.5, ...
        'centerSubregionContourSamples', 40, ...
        'plottedRFoutlineLineWidth', 1.0, ...
        'plottedRFoutlineFaceAlpha', 0.5, ...
        'scaleBarDegs', 1.80);
    % Finalize figure using the Publication-Ready format
    ff.box = 'on';
    PublicationReadyPlotLib.applyFormat(theAxes{1,1},ff);

    % Export figure
    theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    thePDFfileName = fullfile(theRawFiguresDir, 'summaryPlots', 'exampleMRGC.pdf');
    NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);


    [mSequenceStimsNum, nCones] = size(theInputConeMosaicResponses);

    % Compute the temporal equivalent eccentricities in degs
    theTemporalEquivalentEccentricitiesDegs = theMRGCMosaic.temporalEquivalentEccentricityForEccXYDegs(theMRGCMosaic.rgcRFpositionsDegs);

    % Compute the corresponding eccentricities in MMs
    theTemporalEquivalentEccentricitiesMMs = theMRGCMosaic.angularEccDegsToLinearEccMMs(theTemporalEquivalentEccentricitiesDegs);

    fprintf('MRGC mosaic m-sequence RF maps and responses will be saved to %s \n', theMRGCMosaicResponsesFullFileName);
    fprintf('Computing visual m-sequence RF mapping responses for all RGCs in the mosaic ... \n');
    
    
    % Allocate memory
    theMRGCMosaicMSequenceRFmappingResponses = zeros(mSequenceStimsNum, theMRGCMosaic.rgcsNum, 'single');

    % Noise-free responses
    theMRGCMosaic.noiseFlag = 'none';

    % No temporal stuff, just static responses
    nTimePoints = 1; nTrials = 1; 
    theConeMosaicResponseTemporalSupportSeconds = [0];
    theMRGCresponseTemporalSupportSeconds = [];

    parfor iStim = 1:mSequenceStimsNum
        fprintf('Computing mRGC mosaic response for m-sequence frame %d of %d .\n', iStim, mSequenceStimsNum);
        theConeMosaicResponse = reshape(squeeze(theInputConeMosaicResponses(iStim,:)), [nTrials nTimePoints nCones]);
        [theMRGCMosaicResponse, ~, theMRGCresponseTemporalSupportSeconds] = ...
                    theMRGCMosaic.compute(theConeMosaicResponse, theConeMosaicResponseTemporalSupportSeconds);
        theMRGCMosaicMSequenceRFmappingResponses(iStim,:) = single(squeeze(theMRGCMosaicResponse(1,1,:)));
    end % parfor iStim

    % Compute RF maps for all RGCs
    theRFmaps = cell(theMRGCMosaic.rgcsNum, 1);
    parfor iRGC = 1:theMRGCMosaic.rgcsNum
        fprintf('Computing RF for mRGC %d of %d\n', iRGC, theMRGCMosaic.rgcsNum);
        theRFmap = zeros(size(mSequenceIndicatorFunctions,2), size(mSequenceIndicatorFunctions,3));
        allMsequenceStimResponses = theMRGCMosaicMSequenceRFmappingResponses(:, iRGC);
        for iStim = 1:mSequenceStimsNum
            theStimulusFrame = double(squeeze(mSequenceIndicatorFunctions(iStim,:,:)));
            theRFmap = theRFmap + theStimulusFrame * allMsequenceStimResponses(iStim);
        end
        % Expand to full size
        theRFmaps{iRGC} = visualStimulusGenerator.expandFrame(theRFmap, stimParams.rfPixelRetinalPixelsWithin);
    end % iRGC

    % Analyze the maps by fitting a Gaussian ellipsoid
    GaussianEllipsoidSamplesNum = 256;
    theFittedGaussianEllipsoids = RGCMosaicConstructor.helper.simulateExperiment.analyzeSpatialRFmaps(...
        theRFmaps, spatialSupportDegs, GaussianEllipsoidSamplesNum);

    % Update the ellipsoids with sigmas specified in microns
    for iRGC = 1:numel(theFittedGaussianEllipsoids)
        updatedEllipsoid = theFittedGaussianEllipsoids{iRGC};
        % Compute the corresponding sigmas in microns
        updatedEllipsoid.minorSigmaMicrons = 1e3 * theMRGCMosaic.angularEccDegsToLinearEccMMs(updatedEllipsoid.minorSigmaDegs);
        updatedEllipsoid.majorSigmaMicrons = 1e3 * theMRGCMosaic.angularEccDegsToLinearEccMMs(updatedEllipsoid.majorSigmaDegs);
        theFittedGaussianEllipsoids{iRGC} = updatedEllipsoid;
    end


    fprintf('\nSaving computed mRGCRF RF maps and their Gaussian ellipsoid fits to %s ...', theMRGCMosaicResponsesFullFileName);
    save(theMRGCMosaicResponsesFullFileName, ...
            'theMRGCMosaicMSequenceRFmappingResponses', ...
            'theMRGCresponseTemporalSupportSeconds', ...
            'theTemporalEquivalentEccentricitiesDegs', ...
            'theTemporalEquivalentEccentricitiesMMs', ...
            'theRFmaps', ...
            'spatialSupportDegs', ...
            'theFittedGaussianEllipsoids', ...
            '-v7.3');
end


function inputConeMosaicMsequenceRFmapping(theMRGCMosaic, theOI, RFmappingParamsStruct, theInputConeMosaicResponsesFullFileName, visualizeInputConeMosaicResponses)

    fprintf('Computing input cone mosaic M-sequence responses at (x,y) = (%2.1f,%2.1f) using a (%2.2f x %2.2f) deg area\n', ...
        RFmappingParamsStruct.positionDegs(1), RFmappingParamsStruct.positionDegs(2), ...
        RFmappingParamsStruct.sizeDegs(1), RFmappingParamsStruct.sizeDegs(1));


    % Compute cone contrasts for desired chromaticity
    [coneContrasts, totalContrast] = ...
        visualStimulusGenerator.coneContrastsFromChromaticity(RFmappingParamsStruct.chromaticity);

    rfPixelSizeDegs = max(RFmappingParamsStruct.sizeDegs(:))/RFmappingParamsStruct.rfPixelsAcross;
    rfPixelRetinalPixelsWithin = max([1 floor(rfPixelSizeDegs/RFmappingParamsStruct.optimalResolutionDegs)]);
    RFmappingParamsStruct.resolutionDegs = rfPixelSizeDegs/rfPixelRetinalPixelsWithin;
    
    fprintf('\nTo probe RFs using a patch size of %2.3f degs with %d x %d RF pixels\neach RF pixel will be %2.3f degs and contain %d retinal pixels (retinal res of :%2.3f arc min, optimal retinal res: %2.3f arc min).\n', ...
        RFmappingParamsStruct.sizeDegs(1), RFmappingParamsStruct.rfPixelsAcross, RFmappingParamsStruct.rfPixelsAcross, rfPixelSizeDegs, ...
        rfPixelRetinalPixelsWithin, RFmappingParamsStruct.resolutionDegs*60, RFmappingParamsStruct.optimalResolutionDegs*60);

    % Form stimParams struct
    stimParams = struct(...
        'backgroundChromaticity', RFmappingParamsStruct.backgroundChromaticity, ...
        'backgroundLuminanceCdM2', RFmappingParamsStruct.backgroundLuminanceCdM2, ...
        'contrast', totalContrast, ...
        'coneContrasts', coneContrasts, ...
        'sizeDegs', max(RFmappingParamsStruct.sizeDegs(:)), ...
        'positionDegs', RFmappingParamsStruct.positionDegs, ...
        'rfPixelsAcross', RFmappingParamsStruct.rfPixelsAcross, ...
        'rfPixelRetinalPixelsWithin', rfPixelRetinalPixelsWithin, ...
        'rfPixelSizeDegs', rfPixelSizeDegs, ...
        'resolutionDegs', RFmappingParamsStruct.resolutionDegs, ...
        'ternaryInsteadOfBinaryMsequence', RFmappingParamsStruct.ternaryInsteadOfBinaryMsequence, ...
        'mSequenceBitLength', RFmappingParamsStruct.mSequenceBitLength, ...
        'coneMosaicModulationBasedResponse',  true ...
    );

    theMRGCMosaic.visualize(...
        'centerSubregionContourSamples', 40, ...
        'plottedRFoutlineLineWidth', 1.0, ...
        'plottedRFoutlineFaceAlpha', 0.75, ...
        'scaleBarDegs', rfPixelSizeDegs, ...
        'identifyPooledCones', true, ...
        'plotTitle', sprintf('scale bar: %2.3f degs (rf pixel size)', stimParams.rfPixelSizeDegs));

    % Generate presentation display
    viewingDistanceMeters = 4;
    thePresentationDisplay = visualStimulusGenerator.presentationDisplay(...
            theMRGCMosaic.inputConeMosaic.wave, RFmappingParamsStruct.resolutionDegs, ...
            viewingDistanceMeters);

    % Compute the M-sequence spatial patterns
    visualizePatterns = ~true;

    % Compute the m-sequence frame modulation patterns
    mSequenceIndicatorFunctions = visualStimulusGenerator.mSequenceModulationPatterns(...
                    stimParams.rfPixelsAcross, ...
                    'ternaryInsteadOfBinaryMsequence', RFmappingParamsStruct.ternaryInsteadOfBinaryMsequence, ...
                    'mSequenceBitLength', RFmappingParamsStruct.mSequenceBitLength, ...
                    'visualizePatterns', visualizePatterns);
    
    nFrames = size(mSequenceIndicatorFunctions,1);
    theInputConeMosaicResponses = zeros(...
        nFrames, ...
        theMRGCMosaic.inputConeMosaic.conesNum, ...
        'single');
    
    % Compute theNullStimulusScene
    spatialModulationPattern = visualStimulusGenerator.expandFrame(squeeze(mSequenceIndicatorFunctions(1,:,:)), stimParams.rfPixelRetinalPixelsWithin);
    spatialModulationPatternSize = [1 size(spatialModulationPattern,1) size(spatialModulationPattern,2)];
 

    [~, theNullStimulusScene, spatialSupportDegs] = visualStimulusGenerator.stimulusFramesScenes(...
            thePresentationDisplay, stimParams, reshape(spatialModulationPattern, spatialModulationPatternSize), ...
            'frameIndexToCompute', 0, ... % [] field indicates that all stimulus frame scenes must be computed
            'validateScenes', false);

    startingIndex = 1;
    endingIndex = 0;
    while (endingIndex < nFrames)
        if (startingIndex > 1)
            % Only visualize the first batch
            visualizeInputConeMosaicResponses = false;
        end

        endingIndex = min([startingIndex+RFmappingParamsStruct.frameBatchSize-1 nFrames]);
        fprintf('Computing input cone mosaic responses to m-sequence frames %d - %d\n', startingIndex, endingIndex);
        subsequenceFrameIndices = (startingIndex:endingIndex);
        startingIndex = endingIndex+1;

        % Compute theMSequenceFramesSceneSubsequence
        theMSequenceFramesSceneSubsequence = cell(1, numel(subsequenceFrameIndices));

        parfor iFrame = 1:numel(subsequenceFrameIndices)
            theFrameIndex = subsequenceFrameIndices(iFrame);
            % Expand the m-sequence indicator function
            spatialModulationPattern = visualStimulusGenerator.expandFrame(squeeze(mSequenceIndicatorFunctions(theFrameIndex,:,:)), stimParams.rfPixelRetinalPixelsWithin);
            % Compute a scene for this m-sequence frame and add it to theMSequenceFramesSceneSubsequence
            theSceneSequence = visualStimulusGenerator.stimulusFramesScenes(...
                thePresentationDisplay, stimParams, reshape(spatialModulationPattern, spatialModulationPatternSize), ...
                'frameIndexToCompute', 1, ... % [] field indicates that all stimulus frame scenes must be computed
                'validateScenes', false);
            theMSequenceFramesSceneSubsequence{iFrame} = theSceneSequence{1};
        end % iFrame

        % Compute input cone mosaic response to theMSequenceFramesSceneSubsequence
        stimulusPosition = 'mosaic-centered';
        theInputConeMosaicResponses(subsequenceFrameIndices, :) = ...
          RGCMosaicConstructor.helper.simulateExperiment.inputConeMosaicResponseToStimulusFrameSequence(...
            theMRGCMosaic, theOI, theNullStimulusScene, theMSequenceFramesSceneSubsequence, ...
            stimulusPosition, stimParams.coneMosaicModulationBasedResponse, ...
            'visualizeResponse', visualizeInputConeMosaicResponses, ...
            'thePresentationDisplayForVisualizingOpticalSceneOrImage', thePresentationDisplay, ...
            'stimulusInfoString', sprintf('frames: %d - %d', subsequenceFrameIndices(1), subsequenceFrameIndices(end)));
    end % while

    save(theInputConeMosaicResponsesFullFileName, ...
        'theMRGCMosaic', ...
        'RFmappingParamsStruct', ...
        'stimParams', ...
        'mSequenceIndicatorFunctions', ...
        'spatialSupportDegs', ...
        'theInputConeMosaicResponses', ...
        '-v7.3');

end
