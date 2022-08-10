function dryRunRFgeneration()

    % Optics params
    ZernikeDataBase = 'Artal2012';
    examinedSubjectRankOrder = 1;
    pupilDiameterMM = 3.0;

    % Retinal location and eye
    analyzedRetinaMeridian = 'nasal meridian';
    
    % Number of cones in RF center
    conesNumPooledByTheRFcenter = 2;

    analyzedEye = 'right eye';
    subjectRankingEye = 'right eye';

    retinalRFmodelParamsFileName = sprintf('%s_subjRankOrder_%d_%s_%s_pupilDiamMM_%2.2f_conesNumPooledByRFcenter%d.mat', ...
        ZernikeDataBase, examinedSubjectRankOrder, ...
        strrep(analyzedEye, ' ', '_'), ...
        strrep(analyzedRetinaMeridian, ' ', '_'), ...
        pupilDiameterMM, conesNumPooledByTheRFcenter);

    % Load computed data
    load(retinalRFmodelParamsFileName, ...
        'retinalRFparamsDictionary', ...
        'opticsParams', ...
        'targetVisualRFDoGparams');

    % Evaluate generated RFs at a target mosaic eccentricity
    mosaicRadialEccDegs = 0.9;

    % Compute mosaic ecc degs
    [mosaicHorizontalEccDegs, mosaicVerticalEccDegs] = ...
        cMosaic.eccentricitiesForRetinaMeridianInEye(...
            mosaicRadialEccDegs, ...
            analyzedRetinaMeridian, ...
            analyzedEye);
    mosaicEccDegs = [mosaicHorizontalEccDegs mosaicVerticalEccDegs];

    % Retrieve retinal RF params from the saved dictionary at the
    % eccentricity that is closest to the mosaicEcc
    [retinalRFparamsStruct, weightsComputeFunctionHandle, ...
     targetVisualRF, rfSpatialSupportDegs] = dictionaryDataForMosaicEcc(...
        retinalRFparamsDictionary,  opticsParams, mosaicEccDegs);


    mosaicSizeDegs = [1 1];
    theConeMosaic = generateMosaic(analyzedEye, mosaicEccDegs, mosaicSizeDegs);

    targetRFpositionDegs = mosaicEccDegs+0*[0.1 0.3];
    distances = sum((bsxfun(@minus, theConeMosaic.coneRFpositionsDegs, targetRFpositionDegs)).^2,2);
    [~,idx] = sort(distances, 'ascend');

    targetRFCenterConesIndices = idx(1:2);

    if (isequal(weightsComputeFunctionHandle,@RetinaToVisualFieldTransformer.retinalConeWeightsFromDoGmodelParameters))
        pooledConeIndicesAndWeightsStruct = RetinaToVisualFieldTransformer.retinalConeWeightsFromDoGmodelParametersForTargetRFCenterCones(...
            retinalRFparamsStruct, ...
            theConeMosaic, ...
            targetRFCenterConesIndices);
    elseif (isequal(weightsComputeFunctionHandle,@RetinaToVisualFieldTransformer.retinalConeWeightsFromDoGDEmodelParameters))
        error('Implent this')
    else
        error('Uknown function handle');
    end


    % Compute visual RF
    computeVisualRF(pooledConeIndicesAndWeightsStruct, theConeMosaic, opticsParams, targetVisualRF, rfSpatialSupportDegs)

end


function computeVisualRF(pooledConeIndicesAndWeightsStruct, theConeMosaic, opticsParams, targetVisualRF, rfSpatialSupportDegs)
    % Generate optics for the targetConeMosaic
    switch (opticsParams.ZernikeDataBase)
            % Artal
            case RetinaToVisualFieldTransformer.Artal
                subtractCentralRefraction = ArtalOptics.constants.subjectRequiresCentralRefractionCorrection(...
                    theConeMosaic.whichEye, opticsParams.subjID);
                % Polans
            case RetinaToVisualFieldTransformer.Polans
                subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(...
                    theConeMosaic.whichEye, opticsParams.subjID);
    end

    [oiEnsemble, psfEnsemble] = theConeMosaic.oiEnsembleGenerate(...
            theConeMosaic.eccentricityDegs, ...
            'zernikeDataBase', opticsParams.ZernikeDataBase, ...
            'subjectID', opticsParams.subjID, ...
            'pupilDiameterMM', opticsParams.pupilDiameterMM, ...
            'zeroCenterPSF', true, ...
            'subtractCentralRefraction', subtractCentralRefraction, ...
            'wavefrontSpatialSamples', opticsParams.wavefrontSpatialSamples, ...
            'warningInsteadOfErrorForBadZernikeCoeffs', true);


    theOI = oiEnsemble{1};
    thePSFData = psfEnsemble{1};


    % Compute the retinalRF by summing the weighted cone apertures in the
    % center and surround as specified in the computed pooledConeIndicesAndWeightsStruct
    rfSupportX = rfSpatialSupportDegs(:,1);
    rfSupportY = rfSpatialSupportDegs(:,2);
   
    [theRetinalRFcenter, theRetinalRFsurround] = RetinaToVisualFieldTransformer.generateRFsubregionMapsFromPooledCones(...
       rfSupportX, rfSupportY, theConeMosaic, pooledConeIndicesAndWeightsStruct);

    % And the full cone-pooling based retinal RF
    theRetinalRF = theRetinalRFcenter - theRetinalRFsurround;


    % Wavelengths to visualize
    visualizationWavelengths = [450 500 550 600 650];

    figure(1000); clf;
    gain = 1/(abs(max(targetVisualRF(:))));
    
    
    cMap = brewermap(1024, '*RdBu');
    for i = 1:numel(visualizationWavelengths)
        visualizationWavelength = visualizationWavelengths(i);
        [~,iWave] = min(abs(visualizationWavelength-thePSFData.supportWavelength));

        % Select which PSF to use
        theWavePSF = squeeze(thePSFData.data(:,:,iWave));
        
        % Compute the visualRF
        theWaveVisualRF = conv2(theRetinalRF, theWavePSF, 'same');

        ax = subplot(3,numel(visualizationWavelengths), i);
        imagesc(ax, rfSupportX, rfSupportY, gain*theWaveVisualRF/max(abs(targetVisualRF(:))));
        set(ax, 'CLim', [-1 1])
        axis(ax, 'image');
        title(ax, sprintf('%d nm', visualizationWavelength));

        ax = subplot(3,numel(visualizationWavelengths), i + numel(visualizationWavelengths));
        imagesc(ax, rfSupportX, rfSupportY, gain*targetVisualRF/max(abs(targetVisualRF(:))));
        set(ax, 'CLim', [-1 1])
        axis(ax, 'image');
        

        ax = subplot(3,numel(visualizationWavelengths), i + 2*numel(visualizationWavelengths));
        imagesc(ax, rfSupportX, rfSupportY, gain*(theWaveVisualRF-targetVisualRF)/max(abs(targetVisualRF(:))));
        set(ax, 'CLim', [-1 1])
        axis(ax, 'image');
        colormap(cMap);
    end

end


function [retinalRFparamsStruct, weightsComputeFunctionHandle, ...
          targetVisualRF, rfSpatialSupportDegs] = ...
    dictionaryDataForMosaicEcc(retinalRFparamsDictionary, sourceOpticsParams, mosaicEccDegs)

    [sourceHorizontalEccDegs, sourceVerticalEccDegs] = ...
        cMosaic.eccentricitiesForRetinaMeridianInEye(...
            sourceOpticsParams.radialEccDegs, ...
            sourceOpticsParams.analyzedRetinaMeridian, ...
            sourceOpticsParams.analyzedEye);


    % Find the best source data set
    d = sqrt((sourceHorizontalEccDegs-mosaicEccDegs(1)).^2 + ...
             (sourceVerticalEccDegs-mosaicEccDegs(1)).^2);
    [~,iEcc] = min(d(:));
    sourceEccDegs = [sourceHorizontalEccDegs(iEcc) sourceVerticalEccDegs(iEcc)];

    fprintf('Will use the (%2.3f,%2.3f degs) dataset which is closest to the mosaic eccentricity (%2.3f, %2.3f degs)\n', ...
        sourceEccDegs(1), sourceEccDegs(2), mosaicEccDegs(1), mosaicEccDegs(2));

    % Retrieve the dictionary data
    s = retinalRFparamsDictionary(sprintf('Ecc_%2.3f_%2.3f', sourceEccDegs(1), sourceEccDegs(2)));

    retinalRFparamsStruct = s.retinalRFparamsStruct;
    weightsComputeFunctionHandle = s.weightsComputeFunctionHandle;
    targetVisualRF = s.targetVisualRF;
    rfSpatialSupportDegs = s.targetVisualRFspatialSupportDegs;
    maxSpatialSupportDegs = s.maxSpatialSupportDegs;
    theEmployedPSFData = s.theEmployedCircularPSFData;
    clear 's';
end


function theConeMosaic = generateMosaic(analyzedEye, mosaicEccDegs, mosaicSizeDegs)
    theConeMosaic = cMosaic(...
            'whichEye', analyzedEye, ...
            'sizeDegs', mosaicSizeDegs, ...
            'eccentricityDegs', mosaicEccDegs, ...
            'rodIntrusionAdjustedConeAperture', true);
end

