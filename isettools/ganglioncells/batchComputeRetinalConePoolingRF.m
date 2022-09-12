function batchComputeRetinalConePoolingRF(reComputeData)

    % Subjects to use in this batch job
    examinedSubjectRankOrders = 1:38;
    % Remove some subjects which increase the variance a lot
    examinedSubjectRankOrders = setdiff(examinedSubjectRankOrders, [5 16 20 23 31 34 36 37]);


     % Visual RF model to match. Choose between: 
    % {'ellipsoidal gaussian center, gaussian surround', ...
    %  'gaussian center, gaussian surround', ...
    %  'arbitrary center, gaussian surround'}
    visualRFmodel = 'ellipsoidal gaussian center, gaussian surround';


    % Retinal cone pooling model to use. Choose between:
    % 'arbitrary center cone weights, double exponential surround weights'
    % 'arbitrary center cone weights, double gaussian surround weights'
    % 'arbitrary center cone weights, gaussian surround weights'
    % 'arbitrary center cone weights, gaussian surround weights with adjustments'   % takes a long time - not very beneficial 
    retinalConePoolingModel = 'arbitrary center cone weights, double gaussian surround weights';


    for iSubj = 1:numel(examinedSubjectRankOrders)
        examinedSubjectRankOrder = examinedSubjectRankOrders(iSubj);
        batchJob(reComputeData, examinedSubjectRankOrder, visualRFmodel, retinalConePoolingModel);
    end

end

function batchJob(reComputeData, examinedSubjectRankOrder, visualRFmodel, retinalConePoolingModel)
    % Struct with the various optics params
    opticsParams = struct(...
        'rfPositionEccDegs', [0.5 -0.2], ...  % (x,y) eccentricity for the cone pooling RF, in degrees
        'psfCircularSymmetryMode', RetinaToVisualFieldTransformer.psfCircularSymmetryModeAverage, ...
        'ZernikeDataBase', 'Artal2012', ...
        'examinedSubjectRankOrder', examinedSubjectRankOrder, ...
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
        'conesNumPooledByTheRFcenter', 1, ...  % this will depend on the connectivity betwen cone/mRGC mosaics
        'surroundToCenterRcRatio', 6.7, ...
        'surroundToCenterIntegratedSensitivityRatio', 0.55, ...
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
        doDryRunFirst = ~true;
        RetinaToVisualFieldTransformer(theConeMosaic, opticsParams, targetVisualRFDoGparams, ...
            'psfCircularSymmetryMode', opticsParams.psfCircularSymmetryMode, ...
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
