function dryRunRFgeneration()

    % Optics params
    ZernikeDataBase = 'Artal2012';
    examinedSubjectRankOrder = 1;
    pupilDiameterMM = 3.0;

    % Retinal location and eye
    analyzedRetinaMeridian = 'nasal meridian';
    
    % Number of cones in RF center
    conesNumPooledByTheRFcenter = 3;

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

    videoOBJ = VideoWriter(sprintf('ConesNumRFcenter_%d',conesNumPooledByTheRFcenter), 'MPEG-4');
    videoOBJ.FrameRate = 10;
    videoOBJ.Quality = 100;
    videoOBJ.open();

    mosaicRadialEccentricities = [0 0.1 0.25 0.5 0.75 1 1.25 1.5 2 2.5 3 3.5 4 4.5 5:11];
    
    for iEcc = 1:numel(mosaicRadialEccentricities)
        % Evaluate generated RFs at a target mosaic eccentricity
        mosaicRadialEccDegs = mosaicRadialEccentricities(iEcc);

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
         targetVisualRF, rfSpatialSupportDegs, theEmployedPSFData] = dictionaryDataForMosaicEcc(...
            retinalRFparamsDictionary,  opticsParams, mosaicEccDegs);
    
    
        mosaicSizeDegs = mosaicRadialEccDegs*0.4 + 0.4*[1 1];
        theConeMosaic = generateMosaic(analyzedEye, mosaicEccDegs, mosaicSizeDegs);
    
        for iPosition = 1:10
            targetRFpositionDegs = mosaicEccDegs+randn(1,2)*0.15*max(mosaicSizeDegs);
            distances = sum((bsxfun(@minus, theConeMosaic.coneRFpositionsDegs, targetRFpositionDegs)).^2,2);
            [~,idx] = sort(distances, 'ascend');
        
            targetRFCenterConesIndices = idx(1:conesNumPooledByTheRFcenter);
        
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
            theEmployedPSFData = [];
            hFig = computeVisualRF(pooledConeIndicesAndWeightsStruct, theConeMosaic, opticsParams, ...
                theEmployedPSFData, targetVisualRF, rfSpatialSupportDegs);

            drawnow;
            videoOBJ.writeVideo(getframe(hFig));
        end % iPOsition

        %NicePlot.exportFigToPDF(sprintf('analysisEcc%4.3f.pdf', mosaicRadialEccDegs), hFig, 300);
    end
    videoOBJ.close();

        

end


function hFig = computeVisualRF(pooledConeIndicesAndWeightsStruct, theConeMosaic, opticsParams, theEmployedPSFData, targetVisualRF, rfSpatialSupportDegs)
    
    if isempty(theEmployedPSFData)
    
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
    
        % Wavelengths to visualize
        visualizationWavelength = 550;
        [~,iWave] = min(abs(visualizationWavelength-thePSFData.supportWavelength));
        theWavePSF = squeeze(thePSFData.data(:,:,iWave));
    else
        theWavePSF = theEmployedPSFData.data;
    end



    % Compute the retinalRF by summing the weighted cone apertures in the
    % center and surround as specified in the computed pooledConeIndicesAndWeightsStruct
    rfSupportX = rfSpatialSupportDegs(:,1);
    rfSupportY = rfSpatialSupportDegs(:,2);
   
    [theRetinalRFcenter, theRetinalRFsurround] = RetinaToVisualFieldTransformer.generateRFsubregionMapsFromPooledCones(...
       rfSupportX, rfSupportY, theConeMosaic, pooledConeIndicesAndWeightsStruct);

    
    % And the full cone-pooling based retinal RF
    theRetinalRF = theRetinalRFcenter - theRetinalRFsurround;

    
    theVisualRF = conv2(theRetinalRF, theWavePSF, 'same');

    theRetinalRFprofileX = sum(theRetinalRF,1);
    theRetinalRFprofileY = sum(theRetinalRF,2);
    theVisualRFprofileX = sum(theVisualRF,1);
    theVisualRFprofileY = sum(theVisualRF,2);
    theTargetVisualRFprofile = sum(targetVisualRF,1);
    profileRange = [ ...
        min([min(theVisualRFprofileX) min(theVisualRFprofileY) min(theTargetVisualRFprofile)]) ...
        max([max(theVisualRFprofileX) max(theVisualRFprofileY) max(theTargetVisualRFprofile)])];
   
    profileRangeFull = [ ...
        min([min(theVisualRFprofileX) min(theVisualRFprofileY) min(theRetinalRFprofileX) min(theRetinalRFprofileY) min(theTargetVisualRFprofile)]) ...
        max([max(theVisualRFprofileX) max(theVisualRFprofileY) max(theRetinalRFprofileX) max(theRetinalRFprofileY) max(theTargetVisualRFprofile)])];
    deltaP = 0.02*(profileRange(2)-profileRange(1));
    profileRange = profileRange + deltaP*[-1 1];

    idx = find(abs(rfSupportX) < max(rfSupportX)*0.5);
    idy = find(abs(rfSupportY) < max(rfSupportY)*0.5);

    rfSupportX = rfSupportX + theConeMosaic.eccentricityDegs(1);
    rfSupportY = rfSupportY + theConeMosaic.eccentricityDegs(2);

    xRange = rfSupportX(idx);
    xRange = [min(xRange) max(xRange)];
    yRange = rfSupportY(idy);
    yRange = [min(yRange) max(yRange)];

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
        'rowsNum', 2, ...
        'colsNum', 3, ...
        'heightMargin',  0.02, ...
        'widthMargin',    0.03, ...
        'leftMargin',     0.05, ...
        'rightMargin',    0.00, ...
        'bottomMargin',   0.005, ...
        'topMargin',      0.01);

    hFig = figure(1000); clf;
    set(hFig, 'Position', [10 10 1000 1000], 'Color', [1 1 1]);
    ax = subplot('Position', subplotPosVectors(1,1).v);
    plot(ax, rfSupportX, theTargetVisualRFprofile, 'k-', 'LineWidth', 1.5);
    hold on;
    plot(ax, rfSupportX, theTargetVisualRFprofile(:)-0.5*(theVisualRFprofileX(:)+theVisualRFprofileY(:)), 'm--', 'LineWidth', 1.5);
    grid(ax, 'on');
    set(ax, 'YLim', profileRange, 'XLim', xRange, 'FontSize', 16);
    title(ax, 'target visual RF');

    ax = subplot('Position', subplotPosVectors(2,1).v);
    imagesc(ax, rfSupportX, rfSupportY, targetVisualRF);
    axis(ax, 'image');
    grid(ax, 'on');
    set(ax, 'CLim', 1*[-1 1], 'XLim', xRange, 'YLim', yRange, 'YTickLabel', {}, 'FontSize', 16);


    ax = subplot('Position', subplotPosVectors(1,2).v);
    plot(ax, rfSupportX, theTargetVisualRFprofile, 'k-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    hold(ax, 'on')
    plot(ax, rfSupportX, theVisualRFprofileX, 'r-', 'LineWidth', 1.5);
    plot(ax, rfSupportX, theVisualRFprofileY, 'b-', 'LineWidth', 1.5);
    grid(ax, 'on');
    set(ax, 'YLim', profileRange, 'XLim', xRange, 'FontSize', 16);
    title(ax, 'achieved visual RF');
    
    ax = subplot('Position', subplotPosVectors(2,2).v);
    imagesc(ax, rfSupportX, rfSupportY, theVisualRF);
    axis(ax, 'image');
    grid(ax, 'on');
    set(ax, 'CLim', 1*[-1 1], 'XLim', xRange, 'YLim', yRange, 'YTickLabel', {}, 'FontSize', 16);


    ax = subplot('Position', subplotPosVectors(1,3).v);
    plot(ax, rfSupportX, theRetinalRFprofileX, 'r-', 'LineWidth', 1.5);
    hold(ax, 'on')
    plot(ax, rfSupportX, theRetinalRFprofileY, 'b-', 'LineWidth', 1.5);
    grid(ax, 'on');
    set(ax, 'YLim', profileRangeFull, 'XLim', xRange, 'FontSize', 16);
    title(ax, 'retinal RF');

    ax = subplot('Position', subplotPosVectors(2,3).v);
    imagesc(ax, rfSupportX, rfSupportY, theRetinalRF);
    axis(ax, 'image');
    grid(ax, 'on');
    set(ax, 'CLim', 1*[-1 1], 'XLim', xRange, 'YLim', yRange, 'YTickLabel', {}, 'FontSize', 16);

end


function [retinalRFparamsStruct, weightsComputeFunctionHandle, ...
          targetVisualRF, rfSpatialSupportDegs, theEmployedPSFData] = ...
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

