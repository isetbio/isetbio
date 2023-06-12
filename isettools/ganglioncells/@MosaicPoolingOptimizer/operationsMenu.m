function operationSetToPerformContains = operationsMenu(mosaicParams)

    operationDescriptors{1}  = '[ 1] mRGCMosaic generation (step G1)     : Generate the center-connected mRGCMosaic ';
    operationDescriptors{2}  = '[ 2] mRGCMosaic visualization            : Visualize the center-connected mRGCMosaic (and possibly remove unwanted RGCs)';
    operationDescriptors{3}  = '[ 3] mRGCMosaic visualization            : Visualize PSFs at a 3x3 grid within the mRGCMosaic';
    operationDescriptors{4}  = '[ 4] mRGCMosaic generation (step G2)     : Compute STF responses of the input cone mosaic ';
    operationDescriptors{5}  = '[ 5] inputConeMosaic visualization       : Visualize visually-projected cone Rc';

    operationDescriptors{6}  = '[ 6] mRGCMosaic generation (step G3)     : Optimize surround cone pooling models using the input cone mosaic STF responses. (NOTE: This step takes a long time) ';
    operationDescriptors{7}  = '[ 7] mRGCMosaic visualization            : Inspect optimized cone pooling kernels';
    operationDescriptors{8}  = '[ 8] mRGCMosaic generation (step G4)     : Generate the compute-ready mRGCMosaic based on the optimized cone pooling kernels ';
    
    operationDescriptors{9}  = '[ 9] Compute-ready mRGCMosaic validation : Compute visual STFs for all cells in the mosaic';
    operationDescriptors{10} = '[10] Compute-ready mRGCMosaic validation : Fit the DoG model to the computed visual STFs for all cells in the mosaic';
    operationDescriptors{11} = '[11] mRGCMosaic generation (step G5)     : Adjust gain of the compute-ready mRGCMosaic based on the  fitted visual STFs';
    operationDescriptors{12} = '[12] mRGCMosaic generation (step G5)     : Reset gain of the compute-ready mRGCMosaic (1/integrated center weights)';

    operationDescriptors{13} = '[13] Compute-ready mRGCMosaic validation : Visualize cone pooling RF maps and visual STF for individual target RGCs';
    operationDescriptors{14} = '[14] Compute-ready mRGCMosaic validation : Visualize fitted DoG model params for all cells in the mosaic';
    operationDescriptors{15} = '[15] Compute-ready mRGCMosaic validation : Visualize fitted DoG model params for all cells in multiple mosaics';

    operationDescriptors{16} = '[16] Compute-ready mRGCMosaic (RF computation) : Compute visual RFs (subspace) for all cells in the mosaic';
    operationDescriptors{17} = '[17] Compute-ready mRGCMosaic (RF computation) : Visualize visual RF maps (subspace) for individual target RGCs';

    % 01. Generate the center-connected mosaic
    operationSetToPerformContains.generateCenterConnectedRGCMosaic = ~true;

    % 02. Visualize the center-connected mosaic
    operationSetToPerformContains.visualizeCenterConnectedRGCMosaicAndRemoveUnwantedRGCs = ~true;

    % 03. Visualize PSFs at various positions within the mosaic
    operationSetToPerformContains.visualizePSFsWithinRGCMosaic = ~true;

    % 04 Compute cone mosaic STF responses
    operationSetToPerformContains.computeInputConeMosaicSTFresponses = ~true;

    % 05 ComputeInputConeMosaicVisuallyProjectedCharacteristicRadii
    operationSetToPerformContains.computeInputConeMosaicVisuallyProjectedCharacteristicRadii = ~true;

    % 06 Derive optimized surround cone pooling kernels (takes a long time)
    operationSetToPerformContains.optimizeSurroundConePoolingModels = ~true;

    % 07 Inspect optimized cone pooling models
    operationSetToPerformContains.inspectOptimizedSurroundConePoolingModels = ~true;

    % 07 Generate the compute-ready mRGC mosaic
    operationSetToPerformContains.generateComputeReadyMidgetRGCMosaic = ~true;

    % 09 Validation: compute visual STFs for all cells in the mosaic
    operationSetToPerformContains.computeVisualSTFsOfTheComputeReadyMidgetRGCMosaic = ~true;
    
    % 10 Validation: fit the DoG model to the computed visual STFs
    operationSetToPerformContains.fitVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic = ~true;

    % 11 Adjust gain of the compute-ready mRGC mosaic based on fitted visual STFs
    operationSetToPerformContains.adjustGainOfComputeReadyMidgetRGCMosaicBaseOnFittedVisualSTFs = ~true;

    % 12 Reset gain of the compute-ready mRGC mosaic
    operationSetToPerformContains.resetGainOfComputeReadyMidgetRGCMosaic = ~true;

    % 13 Visualization: Cone pooling RF maps and visual STF 
    operationSetToPerformContains.visualizeConePoolingRFmapAndVisualSTFforTargetRGC = ~true;

    % 14 Visualization: Fitted DoG model params for all cells in the current mosaic
    operationSetToPerformContains.visualizeDoGparamsOfVisualSTFsOfSingleMidgetRGCMosaic = ~true;

    % 15 Visualization: Fitted DoG model params for all cells in the current mosaic
    operationSetToPerformContains.visualizeDoGparamsOfVisualSTFsOfMultipleMidgetRGCMosaics = ~true;

    % 16 Subspace RF mapping: compute visual RFs (subspace) for all cells 
    operationSetToPerformContains.computeVisualRFsAcrossTheComputeReadyMidgetRGCMosaic = ~true;

    % 17 Subspace RF mapping: visualize visual RF maps 
    operationSetToPerformContains.visualizeVisualRFmapForTargetRGC = ~true;

    operationSetToPerformContains.animateModelConvergence = ~true;

    invalidActionSelected = true;
    validChoiceIDs = 1:numel(operationDescriptors);

    while (invalidActionSelected)
        % Present options
        fprintf('\n\nAvailable actions for the mosaic at eccentricity (%2.1f,%2.1f):', ...
            mosaicParams.eccDegs(1), mosaicParams.eccDegs(2));
        for iString = 1:numel(operationDescriptors)
            if (strfind(operationDescriptors{iString}, 'step G'))
                fprintf(2,'\n\t%s', operationDescriptors{iString});
            else
                fprintf('\n\t%s', operationDescriptors{iString});
            end

            if (iString == 5) || (iString == 8) || (iString == 12) || (iString == 15)
                fprintf('\n');
            end
        end

        % Get user's choice
        choice = input('\n\nEnter action ID : ', 's');

        if (~isempty(choice))
            choiceID = str2num(choice);
            if (ismember(choiceID, validChoiceIDs))
                switch (choiceID)
                    case 1
                        % Generate center-connected mosaic
                        if (isfield(operationSetToPerformContains, 'generateCenterConnectedRGCMosaic'))
                            operationSetToPerformContains.generateCenterConnectedRGCMosaic = true;
                        else
                             error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'generateCenterConnectedRGCMosaic');
                        end

                    case 2
                        % Visualize center-connected mosaic
                        if (isfield(operationSetToPerformContains, 'visualizeCenterConnectedRGCMosaicAndRemoveUnwantedRGCs'))
                            operationSetToPerformContains.visualizeCenterConnectedRGCMosaicAndRemoveUnwantedRGCs = true;
                        else
                             error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'visualizeCenterConnectedRGCMosaicAndRemoveUnwantedRGCs');
                        end

                    case 3
                        % Visualize PSFs at various positions within the mosaic
                        if (isfield(operationSetToPerformContains, 'visualizePSFsWithinRGCMosaic'))
                            operationSetToPerformContains.visualizePSFsWithinRGCMosaic = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'visualizePSFsWithinRGCMosaic');
                        end

                    case 4
                        % Compute input cone mosaic STF responses
                        if (isfield(operationSetToPerformContains, 'computeInputConeMosaicSTFresponses'))
                            operationSetToPerformContains.computeInputConeMosaicSTFresponses = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'computeInputConeMosaicSTFresponses');
                        end

                    case 5
                        % Compute input cone mosaic visually projected characteristic radii
                        if (isfield(operationSetToPerformContains, 'computeInputConeMosaicVisuallyProjectedCharacteristicRadii'))
                            operationSetToPerformContains.computeInputConeMosaicVisuallyProjectedCharacteristicRadii = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'computeInputConeMosaicVisuallyProjectedCharacteristicRadii');
                        end

                    case 6
                        % Derive optimized surround cone pooling kernels (takes a long time)
                        if (isfield(operationSetToPerformContains, 'optimizeSurroundConePoolingModels'))
                            operationSetToPerformContains.optimizeSurroundConePoolingModels = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'optimizeSurroundConePoolingModels');
                        end

                    case 7
                        % Examine optimized cone pooling models at all or certain  grid nodes
                        if (isfield(operationSetToPerformContains, 'inspectOptimizedSurroundConePoolingModels'))
                            operationSetToPerformContains.inspectOptimizedSurroundConePoolingModels = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'inspectOptimizedSurroundConePoolingModels');
                        end

                    case 8
                        % Generate a compute-ready MRGC mosaic
                        if (isfield(operationSetToPerformContains, 'generateComputeReadyMidgetRGCMosaic'))
                            operationSetToPerformContains.generateComputeReadyMidgetRGCMosaic = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'generateComputeReadyMidgetRGCMosaic');
                        end

                    case 9
                        % Validate a compute-ready mRGCMosaic: step1 - compute visual STFs for all cells
                        if (isfield(operationSetToPerformContains, 'computeVisualSTFsOfTheComputeReadyMidgetRGCMosaic'))
                            operationSetToPerformContains.computeVisualSTFsOfTheComputeReadyMidgetRGCMosaic  = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'computeVisualSTFsOfTheComputeReadyMidgetRGCMosaic');
                        end
              
                    case 10
                        % Validate a compute-ready mRGCMosaic: step2 - fit a DoG model to all the computed visual STFs for all cells
                        if (isfield(operationSetToPerformContains, 'fitVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic'))
                            operationSetToPerformContains.fitVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'fitVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic');
                        end

                    case 11
                        % Adjust gain of compute-read MRGCMosaic based on the computed visual STFs 
                        if (isfield(operationSetToPerformContains, 'adjustGainOfComputeReadyMidgetRGCMosaicBaseOnFittedVisualSTFs'))
                            operationSetToPerformContains.adjustGainOfComputeReadyMidgetRGCMosaicBaseOnFittedVisualSTFs = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'adjustGainOfComputeReadyMidgetRGCMosaicBaseOnFittedVisualSTFs');
                        end

                    case 12
                        % Reset gain of compute-read MRGCMosaic 
                        if (isfield(operationSetToPerformContains, 'resetGainOfComputeReadyMidgetRGCMosaic'))
                            operationSetToPerformContains.resetGainOfComputeReadyMidgetRGCMosaic = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'resetGainOfComputeReadyMidgetRGCMosaic');
                        end

                    case 13
                        % Visualize the derived spatial RF and resulting visual STF for a target RGC 
                        if (isfield(operationSetToPerformContains, 'visualizeConePoolingRFmapAndVisualSTFforTargetRGC'))
                            operationSetToPerformContains.visualizeConePoolingRFmapAndVisualSTFforTargetRGC = true;
                        else
                             error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'visualizeConePoolingRFmapAndVisualSTFforTargetRGC');
                        end

                    case 14
                        % Validate a compute-ready mRGCMosaic: step3 - visualize fitted visual STFs for all cells in a single mRGC mosaic
                        if (isfield(operationSetToPerformContains, 'visualizeDoGparamsOfVisualSTFsOfSingleMidgetRGCMosaic'))
                            operationSetToPerformContains.visualizeDoGparamsOfVisualSTFsOfSingleMidgetRGCMosaic = true;
                        else
                           error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'visualizeDoGparamsOfVisualSTFsOfSingleMidgetRGCMosaic');
                        end

                    case 15
                        % Validate a compute-ready mRGCMosaic: step4 - visualize fitted visual STFs for all cells in multiple mRGC mosaics
                        if (isfield(operationSetToPerformContains, 'visualizeDoGparamsOfVisualSTFsOfMultipleMidgetRGCMosaics'))
                            operationSetToPerformContains.visualizeDoGparamsOfVisualSTFsOfMultipleMidgetRGCMosaics = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'visualizeDoGparamsOfVisualSTFsOfMultipleMidgetRGCMosaics');
                        end
                        
                    case 16
                        if (isfield(operationSetToPerformContains, 'computeVisualRFsAcrossTheComputeReadyMidgetRGCMosaic'))
                            operationSetToPerformContains.computeVisualRFsAcrossTheComputeReadyMidgetRGCMosaic = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'computeVisualRFsAcrossTheComputeReadyMidgetRGCMosaic');
                        end
                        
                    case 17
                        if (isfield(operationSetToPerformContains, 'visualizeVisualRFmapForTargetRGC'))
                            operationSetToPerformContains.visualizeVisualRFmapForTargetRGC  = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: ''%s''.', 'visualizeVisualRFmapForTargetRGC');
                        end
                       
                    otherwise
                        error('Unknown option')
                end % switch

                invalidActionSelected = false;
            else
                fprintf(2, '%d is an invalid option. Choose a number between %d and %d.\n', choiceID, min(validChoiceIDs), max(validChoiceIDs))
            end
        end

    end
   
end
