% Map the visual receptive field of one cone using Hartley patterns
%
% Description:
%    This tutorial illustrates the basic steps of subspace receptive-field
%    (RF) mapping. A compact set of orthogonal Hartley patterns is presented
%    at forward and reverse contrast. The difference between the two cone
%    responses weights each pattern in the reconstructed visual RF.
%
%    For clarity and speed, we map only the cone nearest the mosaic center
%    and use a low-order Hartley basis. The recovered RF is compared with
%    the expected visual RF: the cone aperture convolved with the optical
%    point-spread function (PSF).
%

% History:
%    09/28/22  NPC  ISETBIO Team, Copyright 2022 Wrote it.

function t_conesMapRF

    % Experiment geometry. Eccentricity locates the mosaic on the retina;
    % mosaic size controls how many cones are simulated; stimulus size sets
    % the visual field over which the RF is reconstructed.
    %
    % A small central mosaic is sufficient because this tutorial maps one
    % representative cone rather than producing a movie of every cone.
    mosaicEccDegs = [0 0];
    mosaicSizeDegs = 0.15*[1 1];
    stimSizeDegs = 0.3;

    whichEye = 'right eye';

    % Model the light-collecting aperture of each cone as a Gaussian. These
    % parameters affect the anatomical aperture used in the final comparison.
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

    % Optical parameters. Polans2015 supplies measured wavefronts; the rank
    % chooses a subject from that database; pupil diameter and refractive
    % error determine the PSF. wavefrontSpatialSamples controls PSF sampling
    % accuracy and contributes to startup time, but is computed only once.
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
    [thePSFData, ~, ~, theOI] = computeVlambdaWeightedPSF(...
            opticsParams, theConeMosaic, theConeMosaic.wave);

    % Match the virtual presentation display sampling to the retinal PSF.
    retinalImageResolutionDegs = thePSFData.supportXdegs(2)-thePSFData.supportXdegs(1);
    viewingDistanceMeters = 4;

    theDisplay = rfMappingStimulusGenerator.presentationDisplay(...
        theConeMosaic.wave, retinalImageResolutionDegs, viewingDistanceMeters);

    % Stimulus parameters. coneContrasts=[1 1 1] makes an equal-contrast
    % L/M/S modulation. contrast scales the Hartley pattern around the
    % background. pixelSizeDegs determines the RF-map spatial sampling.
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


    % Hartley spatial patterns used for RF mapping. The number of patterns
    % is (2*omega+1)^2, so omega=1 gives 9 patterns. The previous omega=3
    % produced 49 patterns and required 98 optical-image computations.
    omega = 1;
    HartleySpatialModulationPatterns = ...
        rfMappingStimulusGenerator.HartleyModulationPatterns(...
        omega, stimParams.stimSizeDegs, stimParams.pixelSizeDegs, ...
        'parPoolSize', 1);
    

    % Generate scenes for the Hartley patterns
    [theRFMappingStimulusScenes, theNullStimulusScene, spatialSupportDegs] = rfMappingStimulusGenerator.generateStimulusFramesOnPresentationDisplay(...
        theDisplay, stimParams, HartleySpatialModulationPatterns, ...
        'validateScenes', false);

    fprintf('Generated forward polarity stimuli\n');

    % Generate scenes for the inverse-polarity Hartley patterns
    theReversePolarityRFMappingStimulusScenes = rfMappingStimulusGenerator.generateStimulusFramesOnPresentationDisplay(...
        theDisplay, stimParams, -HartleySpatialModulationPatterns, ...
        'validateScenes', false);
     fprintf('Generated inverse polarity stimuli\n');

    % Select the cone nearest the mosaic center. All cones respond when the
    % mosaic is computed, but only this cone is needed to demonstrate the
    % RF reconstruction.
    coneDistances = vecnorm(theConeMosaic.coneRFpositionsDegs - ...
        mosaicCenterPositionDegs, 2, 2);
    [~, mappedConeIndex] = min(coneDistances);

    % Preallocate the response vector and the RF map for the selected cone.
    nStim = numel(theRFMappingStimulusScenes);
    pixelsNum = size(HartleySpatialModulationPatterns,2);
    mappedConeResponses = zeros(nStim, 1);
    theRFmap = zeros(pixelsNum, pixelsNum);

    % Compute the response to the null scene once, then measure the linear
    % response to each Hartley pattern as forward minus reverse polarity.
    for iFrame = 0:nStim
        fprintf('Computing mosaic response to stim %d of %d\n', iFrame, nStim);

        if (iFrame == 0)
            tic
            theOI = oiCompute(theOI, theNullStimulusScene, 'pad value','mean');
            theNullStimConeMosaicExcitation = squeeze(theConeMosaic.compute(theOI));
            fprintf('Done in %f seconds!!\n', toc);
            continue;
        end

        % Forward-polarity response, expressed as contrast relative to the
        % null response.
        theOI = oiCompute(theOI,theRFMappingStimulusScenes{iFrame},'pad value','mean');
        forwardResponses = (squeeze(theConeMosaic.compute(theOI)) - ...
            theNullStimConeMosaicExcitation)./theNullStimConeMosaicExcitation;
        
        % Reverse-polarity response.
        theOI = oiCompute(theOI,theReversePolarityRFMappingStimulusScenes{iFrame},'pad value','mean');
        reverseResponses = (squeeze(theConeMosaic.compute(theOI)) - ...
            theNullStimConeMosaicExcitation)./theNullStimConeMosaicExcitation;

        % Accumulate the selected cone's RF estimate in the Hartley basis.
        mappedConeResponses(iFrame) = forwardResponses(mappedConeIndex) - ...
            reverseResponses(mappedConeIndex);
        theRFmap = theRFmap + squeeze(...
            HartleySpatialModulationPatterns(iFrame,:,:)) * ...
            mappedConeResponses(iFrame);
    end
    theRFmap = theRFmap/(2*nStim);

    % Summarize the setup and compare the measured RF with its prediction.
    visualizeMappingSetup(theConeMosaic, mappedConeIndex, thePSFData, ...
        spatialSupportDegs, HartleySpatialModulationPatterns, stimParams);
    rfCorrelation = visualizeMappedRF(theConeMosaic, mappedConeIndex, ...
        thePSFData, spatialSupportDegs, theRFmap);

    % A simple quantitative checkpoint for the tutorial.
    fprintf('Mapped cone %d using %d Hartley patterns.\n', ...
        mappedConeIndex, nStim);
    fprintf('Correlation between predicted and mapped RF: %0.3f\n', ...
        rfCorrelation);
    assert(isfinite(rfCorrelation), 'The mapped RF correlation is not finite.');

end

function visualizeMappingSetup(theConeMosaic, mappedConeIndex, thePSFData, ...
    spatialSupportDegs, HartleySpatialModulationPatterns, stimParams)
% Show the mosaic, optics, and one representative mapping stimulus.

    mosaicCenterPositionDegs = mean(theConeMosaic.coneRFpositionsDegs,1);
    xLims = mosaicCenterPositionDegs(1) + round(0.4*stimParams.stimSizeDegs*[-1 1]*100)/100;
    yLims = mosaicCenterPositionDegs(2) + round(0.4*stimParams.stimSizeDegs*[-1 1]*100)/100;

    xTicks = mosaicCenterPositionDegs(1) + round(0.4*stimParams.stimSizeDegs*[-1:0.5:1]*100)/100;
    yTicks = mosaicCenterPositionDegs(2) + round(0.4*stimParams.stimSizeDegs*[-1:0.5:1]*100)/100;


    hFig = figure(1); clf;
    originalFigureUnits = hFig.Units;
    hFig.Units = 'pixels';
    set(hFig, 'Position', [100 300 1200 400], 'Color', [1 1 1]);
    hFig.Units = originalFigureUnits;
    fontSize = 14;

    % Mosaic and the cone whose RF will be mapped.
    axMosaic = subplot(1,3,1);
    theConeMosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', axMosaic, ...
        'domain', 'degrees', ...
        'domainVisualizationLimits', [xLims(1) xLims(2) yLims(1) yLims(2)], ...
        'domainVisualizationTicks', struct('x', xTicks, 'y', yTicks), ...
        'visualizedConeAperture', 'lightCollectingArea4sigma', ...
        'labelConesWithIndices', mappedConeIndex, ...
        'fontSize', fontSize ...
        );
    title(axMosaic, sprintf('mapped cone %d', mappedConeIndex));

    % The wavelength-weighted PSF that blurs the retinal stimulus.
    ax = subplot(1,3,2);
    imagesc(ax, thePSFData.supportXdegs, thePSFData.supportYdegs, thePSFData.data);
    axis(ax, 'image'); axis(ax, 'xy');
    set(ax, 'XLim', xLims-mosaicCenterPositionDegs(1), ...
        'YLim', yLims-mosaicCenterPositionDegs(2), ...
        'FontSize', fontSize);
    colormap(ax, brewermap(1024, 'greys'));
    xlabel(ax, 'space (degrees)');
    ylabel(ax, 'space (degrees)');
    title(ax, 'PSF (vLambda-weighted)');

    % Display a non-DC pattern from the compact Hartley set. The middle
    % pattern is the DC component when omega=1, so use the final pattern.
    representativeStimulusIndex = size(HartleySpatialModulationPatterns,1);
    ax = subplot(1,3,3);
    imagesc(ax, spatialSupportDegs, spatialSupportDegs, ...
        squeeze(HartleySpatialModulationPatterns(representativeStimulusIndex,:,:)));
    axis(ax, 'image'); axis(ax, 'xy');
    set(ax, 'FontSize', fontSize, 'CLim', [-1 1]);
    colormap(ax, brewermap(1024, '*greys'));
    xlabel(ax, 'space (degrees)');
    ylabel(ax, 'space (degrees)');
    title(ax, sprintf('Hartley pattern %d of %d', ...
        representativeStimulusIndex, size(HartleySpatialModulationPatterns,1)));
end


function rfCorrelation = visualizeMappedRF(theConeMosaic, mappedConeIndex, ...
    thePSFData, spatialSupportDegs, theRFmap)
% Compare the anatomical aperture, predicted visual RF, and mapped RF.

    mosaicCenterPositionDegs = mean(theConeMosaic.coneRFpositionsDegs,1);
    conePositionDegs = theConeMosaic.coneRFpositionsDegs(mappedConeIndex,:);
    coneRadiusDegs = theConeMosaic.coneApertureDiametersDegs(mappedConeIndex) * ...
        theConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor;

    % Express the stimulus support relative to the selected cone.
    [xGrid,yGrid] = meshgrid(...
        spatialSupportDegs + mosaicCenterPositionDegs(1) - conePositionDegs(1), ...
        spatialSupportDegs + mosaicCenterPositionDegs(2) - conePositionDegs(2));
    radialDistance = sqrt(xGrid.^2 + yGrid.^2);

    % The anatomical light-collecting aperture and its optically blurred
    % prediction. Normalize each map for shape comparison.
    coneAperture = exp(-(radialDistance/coneRadiusDegs).^2);
    predictedVisualRF = conv2(coneAperture, thePSFData.data, 'same');
    predictedVisualRF = predictedVisualRF/max(predictedVisualRF(:));
    normalizedRFmap = theRFmap/max(abs(theRFmap(:)));

    correlationMatrix = corrcoef(predictedVisualRF(:), normalizedRFmap(:));
    rfCorrelation = correlationMatrix(1,2);

    hFig = figure(2); clf;
    originalFigureUnits = hFig.Units;
    hFig.Units = 'pixels';
    set(hFig, 'Position', [100 100 1200 400], 'Color', [1 1 1]);
    hFig.Units = originalFigureUnits;

    mapTitles = {'anatomical cone aperture', ...
        'predicted visual RF', ...
        sprintf('mapped RF (r = %0.3f)', rfCorrelation)};
    mapData = {coneAperture, predictedVisualRF, normalizedRFmap};

    for ii = 1:numel(mapData)
        ax = subplot(1,3,ii);
        imagesc(ax, spatialSupportDegs, spatialSupportDegs, mapData{ii});
        axis(ax, 'image'); axis(ax, 'xy');
        xlabel(ax, 'space (degrees)');
        ylabel(ax, 'space (degrees)');
        title(ax, mapTitles{ii});
        if ii < 3
            set(ax, 'CLim', [0 1]);
            colormap(ax, brewermap(1024, 'greys'));
        else
            set(ax, 'CLim', [-1 1]);
            colormap(ax, brewermap(1024, '*RdBu'));
        end
    end
end

function [thePSFData, testSubjectID, subtractCentralRefraction, theOI] = computeVlambdaWeightedPSF(...
    opticsParams, theConeMosaic, psfWavelengthSupport)


    % Ensure we have a valid eye specification
    assert(ismember(opticsParams.analyzedEye, {'left eye','right eye'}), ...
        'Invalid analyzed eye specification: ''%s''.', opticsParams.analyzedEye);

    assert(ismember(opticsParams.subjectRankingEye, {'left eye','right eye'}), ...
        'Invalid subject rank eye specification: ''%s''.', opticsParams.subjectRankingEye);


    switch (opticsParams.ZernikeDataBase)
        case 'Artal2012'
            rankedSujectIDs = ArtalOptics.constants.subjectRanking(opticsParams.subjectRankingEye);
            testSubjectID = rankedSujectIDs(opticsParams.examinedSubjectRankOrder);
            subtractCentralRefraction = ArtalOptics.constants.subjectRequiresCentralRefractionCorrection(...
                opticsParams.analyzedEye, testSubjectID);

        case 'Polans2015'
            if (~strcmp(opticsParams.subjectRankingEye, 'right eye'))
                error('Polans measurements exist only for the right eye.');
            end
            rankedSujectIDs = PolansOptics.constants.subjectRanking();
            testSubjectID = rankedSujectIDs(opticsParams.examinedSubjectRankOrder);
            subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(...
                testSubjectID);

        otherwise
            error('Unknown zernike database: ''%ss'.', opticsParams.ZernikeDataBase);
    end


    if (opticsParams.refractiveErrorDiopters == -999)
        % Compute optics at the rf position
        [oiEnsemble, psfEnsemble] = theConeMosaic.oiEnsembleGenerate(opticsParams.positionDegs, ...
                    'zernikeDataBase', opticsParams.ZernikeDataBase, ...
                    'subjectID', testSubjectID, ...
                    'pupilDiameterMM', opticsParams.pupilDiameterMM, ...
                    'refractiveErrorDiopters', 0, ...
                    'zeroCenterPSF', true, ...
                    'subtractCentralRefraction', false, ...
                    'wavefrontSpatialSamples', opticsParams.wavefrontSpatialSamples, ...
                    'upsampleFactor', opticsParams.psfUpsampleFactor, ...
                    'warningInsteadOfErrorForBadZernikeCoeffs', true);
    else
        % Compute optics at the rf position
        [oiEnsemble, psfEnsemble] = theConeMosaic.oiEnsembleGenerate(opticsParams.positionDegs, ...
                    'zernikeDataBase', opticsParams.ZernikeDataBase, ...
                    'subjectID', testSubjectID, ...
                    'pupilDiameterMM', opticsParams.pupilDiameterMM, ...
                    'refractiveErrorDiopters', opticsParams.refractiveErrorDiopters, ...
                    'zeroCenterPSF', true, ...
                    'subtractCentralRefraction', subtractCentralRefraction, ...
                    'wavefrontSpatialSamples', opticsParams.wavefrontSpatialSamples, ...
                    'upsampleFactor', opticsParams.psfUpsampleFactor, ...
                    'warningInsteadOfErrorForBadZernikeCoeffs', true);
    end

    if (isempty(oiEnsemble))
        fprintf(2,'Could not generate optics at this eccentricity');
    end

    % Extract the OTF & the PSF
    thePSFData = psfEnsemble{1};
    theOI = oiEnsemble{1};
    theOTFData.data = theOI.optics.OTF.OTF;
    theOTFData.supportX = theOI.optics.OTF.fx;
    theOTFData.supportY = theOI.optics.OTF.fy;
    theOTFData.supportWavelength = theOI.optics.OTF.wave;

    % Compute v_lambda weights for weigthing the PSF/OTF
    weights = vLambdaWeights(theConeMosaic.wave);

    % Compute vLambda weighted PSF
    vLambdaWeightedPSF = zeros(size(thePSFData.data,1), size(thePSFData.data,2));
    for iWave = 1:size(theOTFData.data,3)
        if (ismember(theConeMosaic.wave(iWave), psfWavelengthSupport)) || ...
           (isempty(psfWavelengthSupport))
            vLambdaWeightedPSF = vLambdaWeightedPSF + thePSFData.data(:,:,iWave) * weights(iWave);
        end
    end
    thePSFData.data = vLambdaWeightedPSF;

    % Ensure we have a unit volume
    thePSFData.data = thePSFData.data / sum(thePSFData.data(:));

    % Specify support in degs instead of the default arc min
    thePSFData.supportXdegs = thePSFData.supportX/60;
    thePSFData.supportYdegs = thePSFData.supportY/60;

    % Remove irrelevant fields
    thePSFData = rmfield(thePSFData, 'supportWavelength');
    thePSFData = rmfield(thePSFData, 'zCoeffs');
    thePSFData = rmfield(thePSFData, 'supportX');
    thePSFData = rmfield(thePSFData, 'supportY');
end

function w = vLambdaWeights(wavelengthSupport)
    load T_xyz1931;
    S = WlsToS(wavelengthSupport(:));
    T_vLambda = SplineCmf(S_xyz1931,T_xyz1931(2,:),S);
    w = T_vLambda/max(T_vLambda(:));
    w = w / sum(w(:));
end
