function [RTVFTobjList, ...
          theOpticsPositionGrid, ...
          theConesNumPooledByTheRFcenterGrid, ...
          theVisualSTFSurroundToCenterRcRatioGrid, ...
          theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid] = generateRTVFobjects(...
                   ZernikeDataBase, subjectRankOrder, pupilDiameterMM, ...
                   theMidgetRGCMosaic, eccentricitySamplingGrid, ...
                   centerConnectableConeTypes, surroundConnectableConeTypes, ...
                   multiStartsNum)

    % Find out the range of cones in the RF center
    allConesNumPooledByTheRFcenters = full(sum(theMidgetRGCMosaic.rgcRFcenterConeConnectivityMatrix,1));
    conesNumPooledByTheRFcenters = unique(allConesNumPooledByTheRFcenters);
    fprintf('Cones/RF center for this mosaic: %d\n', conesNumPooledByTheRFcenters);

    % Generate grid of optical positions x conesNumInRFcenter
    [RTVTOpticalPositionGrid, RTVTConesNumInRFcenterGrid] = ...
        meshgrid(1:size(eccentricitySamplingGrid,1), conesNumPooledByTheRFcenters);
    

    % Allocate memory
    RTVTobjectsNum = numel(RTVTConesNumInRFcenterGrid);
    RTVFTobjList = cell(RTVTobjectsNum,1);
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

        % Get indices of center cones for the RGC that is closest to the opticsPositionDegs
        indicesOfRGCsWithThisManyCenterCones  = find(allConesNumPooledByTheRFcenters == conesNumPooled);
        [~,idx] = min(sum((bsxfun(@minus, theMidgetRGCMosaic.rgcRFpositionsDegs(indicesOfRGCsWithThisManyCenterCones,:), opticsPositionDegs)).^2,2));
        theTargetRGCindex = indicesOfRGCsWithThisManyCenterCones(idx);
        indicesOfConesPooledByTheRFcenter = find(theMidgetRGCMosaic.rgcRFcenterConeConnectivityMatrix(:,theTargetRGCindex)> 0);
        typesOfConesPooledByTheRFcenter = theMidgetRGCMosaic.inputConeMosaic.coneTypes(indicesOfConesPooledByTheRFcenter);

        % Assert that we have the correct number of center cones
        assert((numel(indicesOfConesPooledByTheRFcenter) == conesNumPooled), ...
            sprintf('indicesOfConesPooledByTheRFcenter should have %d entries but it has %d.', numel(indicesOfConesPooledByTheRFcenter), conesNumPooled));
        
        % Assert that the center cones are all connectable
        assert(all(ismember(typesOfConesPooledByTheRFcenter, centerConnectableConeTypes)), ...
            sprintf('indicesOfConesPooledByTheRFcenter are not all connectable'));

        % Report types of cones in the RF center
        if (conesNumPooled == 1)
            fprintf('\nCone type in single-cone RF center: ');
        else
            fprintf('\nCone types in multi-cone RF center: ');
        end

        for iInputConeIndex = 1:numel(typesOfConesPooledByTheRFcenter)
            switch (typesOfConesPooledByTheRFcenter(iInputConeIndex))
                case cMosaic.LCONE_ID
                    fprintf('L ');
                case cMosaic.MCONE_ID
                    fprintf('M ');
                case cMosaic.SCONE_ID
                    fprintf('S ');
                otherwise
                    error('not an LMS cone type');
            end
        end
        fprintf('\n');

        % Unit weights for all center cones
        weightsOfConesPooledByTheRFcenter = ones(1,numel(indicesOfConesPooledByTheRFcenter));

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
            'centerConnectableConeTypes', centerConnectableConeTypes, ...
            'surroundConnectableConeTypes', surroundConnectableConeTypes, ...
            'indicesOfConesPooledByTheRFcenter', indicesOfConesPooledByTheRFcenter, ...
            'weightsOfConesPooledByTheRFcenter', weightsOfConesPooledByTheRFcenter, ...
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
