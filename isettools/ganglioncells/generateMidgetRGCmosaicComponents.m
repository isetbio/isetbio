function generateMidgetRGCmosaicComponents
    eccDegs = [0 0];
    sizeDegs = [0.5 0.5];

    fName = sprintf('mRGCmosaicComponents_eccDegs_%2.2f.mat', eccDegs(1));

    % Operations
    operations = {...
        'generateMosaic' ...
        'generateRTVFTobjList' ...
        'wireRFsurrounds' ...
        };

    % Generate the RTVFTobjList
    operations = operations(2);

    for iOp = 1:numel(operations)

        switch (operations{iOp})
            case 'generateMosaic'
                fprintf('Generating midgetRGCmosaic ...\n');
                % Generate the midgetRGCmosaic
                theMidgetRGCmosaic = midgetRGCMosaic(...
                        'sourceLatticeSizeDegs', 60, ...
                        'eccentricityDegs', eccDegs, ...
                        'sizeDegs', sizeDegs ...
                        );
                save(fName, 'theMidgetRGCmosaic', '-v7.3');
                fprintf('Exported the computed midgetRGCMosaic to %s\n', fName);

            case 'generateRTVFTobjList'
                fprintf('Generating RTVFTobjList ... \n');

                % Load the midget RGCmosaic
                load(fName, 'theMidgetRGCmosaic');

                % Generate the list of RTVFT objects for different gridPositions
                iGridPosition = 1;

                eccDegsGrid(iGridPosition,:) = [0 0];
                conesNumPooledByTheRFcenterGrid(iGridPosition) = 1;

                % From Croner & Kaplan '95 (Figure 4c and text)
                % "P surrounds were on average 6.7 times wider than the centers of
                % the same cells, or about 45 times larger in area".
                surroundToCenterRcRatioGrid(iGridPosition) = 6.7;

                % Also from Croner & Kaplan '95 (Figure 10b)
                % "These mean ratios for P and M cells are not significantly different
                % (Student's t-test: P = 0.482). The overall mean ratio is 0.55.
                surroundToCenterIntegratedSensitivityRatioGrid(iGridPosition) = 0.55;
                

                % Compute a list of RTVFTobj for each of the examined grid positions
                RTVFTobjList = generateRTVFTobjects(theMidgetRGCmosaic, ...
                    eccDegsGrid, conesNumPooledByTheRFcenterGrid, ...
                    surroundToCenterRcRatioGrid, surroundToCenterIntegratedSensitivityRatioGrid);

                % Save the computed list of RTVFTobj for each of the examined grid positions
                save(fName, ...
                    'RTVFTobjList', ...
                    'eccDegsGrid', ...
                    'conesNumPooledByTheRFcenterGrid', ...
                    'surroundToCenterRcRatioGrid', ...
                    'surroundToCenterIntegratedSensitivityRatioGrid', ...
                    '-append');
                fprintf('Appended the computed RTVFTobjList to %s\n', fName);

            case 'wireRFsurrounds'
                load(fName, ...
                    'theMidgetRGCmosaic', ...
                    'RTVFTobjList', ...
                    'eccDegsGrid', ...
                    'conesNumPooledByTheRFcenterGrid', ...
                    'surroundToCenterRcRatioGrid', ...
                    'surroundToCenterIntegratedSensitivityRatioGrid');

                    RTVFTobjList{1}

            otherwise
                error('Unknown operation: ''%s''.', operations{iOp});
        end % Switch
    end
end

function RTVFTobjList = generateRTVFTobjects(theMidgetRGCmosaic, ...
    eccDegsGrid, conesNumPooledByTheRFcenterGrid, ...
    surroundToCenterRcRatioGrid, surroundToCenterIntegratedSensitivityRatioGrid)

    gridPositionsNum = size(eccDegsGrid,1);
    RTVFTobjList = cell(1, gridPositionsNum);

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


    targetVisualRFDoGparams = struct(...
        'conesNumPooledByTheRFcenter', [], ...  % this will depend on the connectivity betwen cone/mRGC mosaics
        'surroundToCenterRcRatio', [], ...
        'surroundToCenterIntegratedSensitivityRatio', [], ... 
        'visualRFmodel', visualRFmodel, ...  
        'retinalConePoolingModel', retinalConePoolingModel ...
        );

    ZernikeDataBase = 'Artal2012';
    subjectRankOrder = 2;

    % Struct with the various optics params
    opticsParams = struct(...
        'positionDegs', [], ...           % (x,y) eccentricity for the PSF, in degrees
        'ZernikeDataBase', ZernikeDataBase, ...
        'examinedSubjectRankOrder', subjectRankOrder, ...
        'refractiveErrorDiopters', 0.0, ...    % use -999 for optics that do not subtract the central refraction
        'analyzedEye', 'right eye', ...
        'subjectRankingEye', 'right eye', ...
        'pupilDiameterMM', 3.0, ...
        'wavefrontSpatialSamples', 401, ...
        'psfUpsampleFactor', 2 ...
        );

    % Extract the cone mosaic from the midgetRGCmosaic
    theConeMosaic = theMidgetRGCmosaic.inputConeMosaic;

    parfor iGridPosition = 1:gridPositionsNum
        % Copy params structs
        theGridOpticsParams = opticsParams;
        theGridTargetVisualRFDoGparams = targetVisualRFDoGparams;
        
        % Update params structs for this grid position

        % Update opticsParams position for this grid position
        theGridOpticsParams.positionDegs = eccDegsGrid(iGridPosition,:);

        % Update targetVisualRFDoGparams conesNum for this grid position
        theGridTargetVisualRFDoGparams.conesNumPooledByTheRFcenter = conesNumPooledByTheRFcenterGrid(iGridPosition);
        
        % Update targetVisualRFDoGparams surroundToCenterRcRatio for this grid position
        theGridTargetVisualRFDoGparams.surroundToCenterRcRatio = surroundToCenterRcRatioGrid(iGridPosition);

        % Update targetVisualRFDoGparams surroundToCenterIntegratedSensitivityRatio
        theGridTargetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio = surroundToCenterIntegratedSensitivityRatioGrid(iGridPosition);

        % Compute the RetinaToVisualFieldTransformer for this grid position
        tic
        multiStartsNum = 2;
        doDryRunFirst = true;
        RTVFTobjList{iGridPosition} = RetinaToVisualFieldTransformer(...
            theConeMosaic, ...
            theGridOpticsParams, theGridTargetVisualRFDoGparams, ...
            'simulateCronerKaplanEstimation', true, ...
            'multiStartsNum', multiStartsNum, ...
            'doDryRunFirst', doDryRunFirst);
        toc

    end % iGridPosition

end
