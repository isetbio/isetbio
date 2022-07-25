function demoVisualFieldToRetinaDeconvolution

    xFormer = RetinaToVisualFieldTransformer('ZernikeDataBase', 'Artal2012');

    analyzedRetinaMeridian = 'nasal meridian';
    subjectRankingEye = 'left eye';
    analyzedEye = 'left eye';
    pupilDiameterMM = 3.0;

    radialEccDegs = 3; %0:1:30;


    % Only the first 30 subjects
    examinedSubjectRankOrders = 1:38;
    % Remove some subjects which increase the
    examinedSubjectRankOrders = setdiff(examinedSubjectRankOrders, [5 16 20 23 31 34 36 37]);
    examinedSubjectRankOrders = 2;

    % Get the horizontal and vertical eccs corresponding to the target radial ecc at the
    % analyzed retinal meridian and eye
    
    [horizontalEccDegs, verticalEccDegs] = cMosaic.eccentricitiesForRetinaMeridianInEye(...
            radialEccDegs, analyzedRetinaMeridian, analyzedEye);

    % Struct with DoG model params for a midget RGC
    visualRFDoGparams = struct(...
        'Kc', 1, ...
        'RcDegs', [], ...  % empty indicates that we should derive this based on the cone mosaic and conesNumInRFcenter
        'surroundToCenterRcRatio', 3, ...
        'surroundToCenterIntegratedRatio', 0.8);

    % How many cones are pooled by 
    % the center mechanism
    conesNumPooledByTheRFcenter = 2;

    % Reduce high-frequency deconvolution noise
    regularizationAlphaRange = [0.001 0.5];

    % Analyze all examined subjects
    for iSubj = 1:numel(examinedSubjectRankOrders)
        % Get the subject ID
        subjectRankOrder = examinedSubjectRankOrders(iSubj);
        subjID = xFormer.subjectWithRankInEye(subjectRankOrder, subjectRankingEye);

        for iEcc = 1:numel(horizontalEccDegs)
            % Analyze effect of optics at this eccentricity
            eccDegs = [horizontalEccDegs(iEcc) verticalEccDegs(iEcc)];
            
            deconvolvedVisualRFdata = xFormer.deconvolveVisualRFgivenOptics(...
                visualRFDoGparams, conesNumPooledByTheRFcenter, eccDegs, analyzedEye, subjID, pupilDiameterMM, ...
                regularizationAlphaRange);
        end % iEcc
    end % iSubj

end