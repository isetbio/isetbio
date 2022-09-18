function batchComputeRetinalConePoolingRF(reComputeData)

    % Subjects to use in this batch job
    examinedSubjectRankOrders = 1:41;
    % Remove some subjects which increase the variance a lot
    examinedSubjectRankOrders = setdiff(examinedSubjectRankOrders, [5 16 20 23 31 34 36 37]);


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
 %   retinalConePoolingModel = 'arbitrary center cone weights, gaussian surround weights';

%     for iSubj = 1:numel(examinedSubjectRankOrders)
%         examinedSubjectRankOrder = examinedSubjectRankOrders(iSubj);
%         batchJob(reComputeData, examinedSubjectRankOrder, visualRFmodel, retinalConePoolingModel);
%     end


%      examinedHorizontalEccDegs = [2];
%      conesNumPooledByTheRFcenter = 1;
%      examinedSubjectRankOrder = 2;
% 
%     for iEcc = 1:numel(examinedHorizontalEccDegs)
%         batchJob(reComputeData, examinedSubjectRankOrder, examinedHorizontalEccDegs(iEcc), conesNumPooledByTheRFcenter, visualRFmodel, retinalConePoolingModel);
%     end

%      examinedHorizontalEccDegs = 5;
%      conesNumPooledByTheRFcenterList = [2 3 4 5];
%      examinedSubjectRankOrder = 2;
% % 
%      for iCones = 1:numel(conesNumPooledByTheRFcenterList)
%          batchJob(reComputeData, examinedSubjectRankOrder, examinedHorizontalEccDegs, conesNumPooledByTheRFcenterList(iCones), visualRFmodel, retinalConePoolingModel);
%      end

    comparison = 'acrossSubjects';
    comparison = 'acrossEccentricities';
    comparison = 'acrossNumberOfConesInRFcenter';
    comparison = 'acrossRFspatialCharacteristics & surround cone pooling models';


    switch (comparison)
        case 'acrossSubjects'
            examinedHorizontalEccDegs = 0.5;
            conesNumPooledByTheRFcenter = 1;
            
            surroundToCenterRcRatio = 6.7;
            surroundToCenterIntegratedSensitivityRatio = 0.55;

            for iSubj = 1:4:numel(examinedSubjectRankOrders)
                examinedSubjectRankOrder = examinedSubjectRankOrders(iSubj);
                batchJob(reComputeData, examinedSubjectRankOrder, examinedHorizontalEccDegs, ...
                    surroundToCenterRcRatio, surroundToCenterIntegratedSensitivityRatio, ...
                    conesNumPooledByTheRFcenter, visualRFmodel, retinalConePoolingModel);
            end

        case 'acrossEccentricities'
            examinedHorizontalEccDegsList = [0.5 1 2 4 8];
            examinedSubjectRankOrder = 41;
            conesNumPooledByTheRFcenter = 1;

            surroundToCenterRcRatio = 6.7;
            surroundToCenterIntegratedSensitivityRatio = 0.55;

            for iEcc = 1:numel(examinedHorizontalEccDegsList)
                examinedHorizontalEccDegs = examinedHorizontalEccDegsList(iEcc);
                batchJob(reComputeData, examinedSubjectRankOrder, examinedHorizontalEccDegs, ...
                    surroundToCenterRcRatio, surroundToCenterIntegratedSensitivityRatio, ...
                    conesNumPooledByTheRFcenter, visualRFmodel, retinalConePoolingModel);
            end

        case 'acrossNumberOfConesInRFcenter'
            examinedHorizontalEccDegs = 2;
            examinedSubjectRankOrder = 29;

            surroundToCenterRcRatio = 6.7;
            surroundToCenterIntegratedSensitivityRatio = 0.55;

            conesNumPooledByTheRFcenterList = [1 2 3 4 5 6 7];
            for iConesNum = 1:numel(conesNumPooledByTheRFcenterList)
                conesNumPooledByTheRFcenter = conesNumPooledByTheRFcenterList(iConesNum);
                batchJob(reComputeData, examinedSubjectRankOrder, examinedHorizontalEccDegs, ...
                    surroundToCenterRcRatio, surroundToCenterIntegratedSensitivityRatio, ...
                    conesNumPooledByTheRFcenter, visualRFmodel, retinalConePoolingModel);
            end

        case 'acrossRFspatialCharacteristics & surround cone pooling models'
            examinedCSsensitivityRatiosList = [0.55 0.2 0.2 0.8 0.8];
            examinedRsRcRatiosList          = [6.7    3 10    3 10];
            conesNumPooledByTheRFcenter = 1;

            examinedHorizontalEccDegs = 0.25;
            examinedSubjectRankOrder = 2;
            retinalConePoolingModel = 'arbitrary center cone weights, gaussian surround weights';

            for iShape = 1:numel(examinedCSsensitivityRatiosList)
                surroundToCenterIntegratedSensitivityRatio = examinedCSsensitivityRatiosList(iShape);
                surroundToCenterRcRatio = examinedRsRcRatiosList(iShape);
                batchJob(reComputeData, examinedSubjectRankOrder, examinedHorizontalEccDegs, ...
                    surroundToCenterRcRatio, surroundToCenterIntegratedSensitivityRatio, ...
                    conesNumPooledByTheRFcenter, visualRFmodel, retinalConePoolingModel);
            end

            retinalConePoolingModel = 'arbitrary center cone weights, double exponential surround weights-free';

            for iShape = 1:numel(examinedCSsensitivityRatiosList)
                surroundToCenterIntegratedSensitivityRatio = examinedCSsensitivityRatiosList(iShape);
                surroundToCenterRcRatio = examinedRsRcRatiosList(iShape);
                batchJob(reComputeData, examinedSubjectRankOrder, examinedHorizontalEccDegs, ...
                    surroundToCenterRcRatio, surroundToCenterIntegratedSensitivityRatio, ...
                    conesNumPooledByTheRFcenter, visualRFmodel, retinalConePoolingModel);
            end

    end % switch


end


function batchJob(reComputeData, examinedSubjectRankOrder, horizontalEccDegs, ...
    surroundToCenterRcRatio, surroundToCenterIntegratedSensitivityRatio, conesNumPooledByTheRFcenter, ...
    visualRFmodel, retinalConePoolingModel)

    % Struct with the various optics params
    opticsParams = struct(...
        'rfPositionEccDegs', [horizontalEccDegs 0.0], ...  % (x,y) eccentricity for the cone pooling RF, in degrees
        'ZernikeDataBase', 'Artal2012', ...
        'examinedSubjectRankOrder', examinedSubjectRankOrder, ...
        'refractiveErrorDiopters', 0.0, ...   % use -999 for optics that do not subtract the central refraction
        'analyzedEye', 'right eye', ...
        'subjectRankingEye', 'right eye', ...
        'pupilDiameterMM', 3.0, ...
        'wavefrontSpatialSamples', 801 ...
        );
  
    % From Croner & Kaplan '95 (Figure 4c and text)
    % "P surrounds were on average 6.7 times wider than the centers of
    % the same cells, or about 45 times larger in area".
    
    % Also from Croner & Kaplan '95 (Figure 10b)
    % "These mean ratios for P and M cells are not significantly different
    % (Student's t-test: P = 0.482). The overall mean ratio is 0.55.

    targetVisualRFDoGparams = struct(...
        'conesNumPooledByTheRFcenter', conesNumPooledByTheRFcenter, ...  % this will depend on the connectivity betwen cone/mRGC mosaics
        'surroundToCenterRcRatio', surroundToCenterRcRatio, ...
        'surroundToCenterIntegratedSensitivityRatio', surroundToCenterIntegratedSensitivityRatio, ... 
        'visualRFmodel', visualRFmodel, ...  
        'retinalConePoolingModel', retinalConePoolingModel ...
        );
    
    if (reComputeData)

        % The input cone mosaic 
        % Use a Gaussian cone aperture with
        % sigma equal to 0.204 x inner segment diameter (cone diameter)
        sigmaGaussian = 0.204;  % From McMahon et al, 2000
        
        % Set cone aperture modifiers
        coneApertureModifiers = struct(...
            'smoothLocalVariations', true, ...
            'sigma',  sigmaGaussian, ...
            'shape', 'Gaussian');
        
        theConeMosaic = cMosaic(...
            'sourceLatticeSizeDegs', 60, ...
            'whichEye', opticsParams.analyzedEye, ...
            'eccentricityDegs', opticsParams.rfPositionEccDegs, ...
            'sizeDegs', [1 1], ...
            'rodIntrusionAdjustedConeAperture', true, ...
            'coneApertureModifiers', coneApertureModifiers, ...
            'customDegsToMMsConversionFunction', @RGCmodels.Watson.convert.rhoDegsToMMs, ...
            'customMMsToDegsConversionFunction', @RGCmodels.Watson.convert.rhoMMsToDegs);


        % Instantiate the RetinaToVisualFieldTransformer
        tic
        multiStartsNum = 16;
        doDryRunFirst = true;
        RetinaToVisualFieldTransformer(theConeMosaic, opticsParams, targetVisualRFDoGparams, ...
            'simulateCronerKaplanEstimation', true, ...
            'multiStartsNum', multiStartsNum, ...
            'doDryRunFirst', doDryRunFirst);
        toc

    else
        % Load a previously saved computed object and visualize the results
        computedObjectDataFileName = RetinaToVisualFieldTransformer.computedObjectDataFileName(opticsParams, targetVisualRFDoGparams);
        d = load(computedObjectDataFileName);
        obj = d.obj; clear 'd';

        % Visualize the results
        obj.visualizeResults();

        % Visualize the fitted params
        RetinaToVisualFieldTransformer.visualizeFittedParamValues(obj.rfComputeStruct.retinalConePoolingParams);
    end


end
