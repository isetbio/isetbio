function [RTVFTobjList, ...
          theOpticsPositionGrid, ...
          theConesNumPooledByTheRFcenterGrid, ...
          theVisualSTFSurroundToCenterRcRatioGrid, ...
          theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid] = generateRTVFobjects(...
                   ZernikeDataBase, subjectRankOrder, pupilDiameterMM, ...
                   theMidgetRGCMosaic, eccentricitySamplingGrid, ...
                   multiStartsNum)

    % Find out the range of cones in the RF center
    conesNumPooledByTheRFcenters = unique(full(sum(theMidgetRGCMosaic.rgcRFcenterConeConnectivityMatrix,1)));
    fprintf('Cones/RF center for this mosaic: %d\n', conesNumPooledByTheRFcenters);

    % Generate grid of optical positions x conesNumInRFcenter
    [RTVTOpticalPositionGrid, RTVTConesNumInRFcenterGrid] = ...
        meshgrid(1:size(eccentricitySamplingGrid,1), conesNumPooledByTheRFcenters);
    

    % Allocate memory
    RTVTobjectsNum = numel(RTVTConesNumInRFcenterGrid);
    RTVFTobjList = cell(1, RTVTobjectsNum);
    theOpticsPositionGrid = zeros(RTVTobjectsNum,2);
    theConesNumPooledByTheRFcenterGrid = zeros(RTVTobjectsNum,1);
    theVisualSTFSurroundToCenterRcRatioGrid = zeros(RTVTobjectsNum,1);
    theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid = zeros(RTVTobjectsNum,1);

    for iRTVobjIndex = 1:RTVTobjectsNum
        % Optics position
        eccPositionIndex = RTVTOpticalPositionGrid(iRTVobjIndex);
        opticsPositionDegs = eccentricitySamplingGrid(eccPositionIndex,:);

        % Cones num in RF center
        conesNumPooled = RTVTConesNumInRFcenterGrid(iRTVobjIndex);

        % Rs/Rc ratio
        surroundToCenterRcRatio = RGCmodels.CronerKaplan.constants.surroundToCenterRcRatio;

        % Temporal-equivalent eccentricity based SCint sensitivity ratio
        temporalEquivalentEccDegs = theMidgetRGCMosaic.temporalEquivalentEccentricityForEccentricity(opticsPositionDegs);
        radialTemporalEquivalentEccDegs = sqrt(sum(temporalEquivalentEccDegs.^2,2));
        scIntSensitivity = RGCmodels.CronerKaplan.constants.surroundToCenterIntegratedSensitivityRatioFromEccDegsForPcells(radialTemporalEquivalentEccDegs);

        % Optics params
        [psfUpsampleFactor, wavefrontSpatialSamples] = psfAndWavefrontSamples(radialTemporalEquivalentEccDegs);
        theOpticsParams = struct(...
            'ZernikeDataBase', ZernikeDataBase, ...
            'examinedSubjectRankOrder', subjectRankOrder, ...
            'pupilDiameterMM', pupilDiameterMM, ...
            'positionDegs', opticsPositionDegs, ...
            'refractiveErrorDiopters', 0.0, ...    % use -999 for optics that do not subtract the central refraction
            'analyzedEye', theMidgetRGCMosaic.whichEye, ...
            'subjectRankingEye', 'right eye', ...
            'psfUpsampleFactor', psfUpsampleFactor, ...
            'wavefrontSpatialSamples', wavefrontSpatialSamples ...
        );

        % TargetVisualRFDoGparams

        % Visual RF model to match. Choose between: 
        % {'ellipsoidal gaussian center, gaussian surround', ...
        %  'gaussian center, gaussian surround', ...
        %  'arbitrary center, gaussian surround'}
        % When simulating the Croner&Kaplan assessment this must be set to 'gaussian center, gaussian surround';
        visualRFmodel = 'gaussian center, gaussian surround';
    
        % Retinal cone pooling model to use. Choose between:
        % 'arbitrary center cone weights, variable exponential surround weights';
        % 'arbitrary center cone weights, double exponential surround weights-free'
        % 'arbitrary center cone weights, double exponential surround weights-meanVnVwRatio'
        % 'arbitrary center cone weights, double exponential surround weights-meanRnRwRatio'
        % 'arbitrary center cone weights, double gaussian surround weights'
        % 'arbitrary center cone weights, gaussian surround weights'
        % 'arbitrary center cone weights, gaussian surround weights with adjustments'   % takes a long time - not very beneficial 
        % retinalConePoolingModel = 'arbitrary center cone weights, double exponential surround weights-meanRnRwRatio';
        retinalConePoolingModel = 'arbitrary center cone weights, double exponential surround weights-free';

        theTargetVisualRFDoGparams = struct(...
            'visualRFmodel', visualRFmodel, ... 
            'retinalConePoolingModel', retinalConePoolingModel, ...
            'conesNumPooledByTheRFcenter', conesNumPooled, ...
            'surroundToCenterRcRatio', surroundToCenterRcRatio, ...
            'surroundToCenterIntegratedSensitivityRatio', scIntSensitivity);

        

        % Compute the RetinaToVisualFieldTransformer for this (optics position, conesNumPooled) pair
        RTVFTobjList{iRTVobjIndex} = RetinaToVisualFieldTransformer(...
            theMidgetRGCMosaic.inputConeMosaic, ...
            theOpticsParams, ...
            theTargetVisualRFDoGparams, ...
            'simulateCronerKaplanEstimation', true, ...
            'multiStartsNum', multiStartsNum , ...
            'doDryRunFirst', true);


        theOpticsPositionGrid(iRTVobjIndex,:) = opticsPositionDegs;
        theConesNumPooledByTheRFcenterGrid(iRTVobjIndex) = conesNumPooled;
        theVisualSTFSurroundToCenterRcRatioGrid(iRTVobjIndex) = surroundToCenterRcRatio;
        theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid(iRTVobjIndex) = scIntSensitivity;
    end
end

function [psfUpsampleFactor, wavefrontSpatialSamples] = psfAndWavefrontSamples(radialEccDegs)
    if (radialEccDegs <= 0.5)
        % Cones are tiny, so upsample the spatial resolution
        psfUpsampleFactor = 2;
        wavefrontSpatialSamples = 301;
    elseif (radialEccDegs <= 8)
        psfUpsampleFactor = 1;
        wavefrontSpatialSamples = 401;
    elseif (radialEccDegs <= 14)
        psfUpsampleFactor = 1;
        wavefrontSpatialSamples = 601;   
    else
        psfUpsampleFactor = 1;
        wavefrontSpatialSamples = 701;
    end

end
