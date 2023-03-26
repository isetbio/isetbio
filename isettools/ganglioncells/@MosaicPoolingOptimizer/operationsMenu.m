function operationSetToPerformContains = operationsMenu(mosaicParams)

    operationDescriptors{1}  = '[ 1] mRGCMosaic generation (step 1) : Generate the center-connected mRGCMosaic';
    operationDescriptors{2}  = '[ 2] mRGCMosaic visualization       : Visualize the center-connected mRGCMosaic (and possibly remove unwanted RGCs)';
    operationDescriptors{3}  = '[ 3] mRGCMosaic visualization       : Visualize PSFs at a 3x3 grid within the mRGCMosaic';
    operationDescriptors{4}  = '[ 4] mRGCMosaic generation (step 2) : Compute STF responses of the input cone mosaic';
    operationDescriptors{5}  = '[ 5] mRGCMosaic generation (step 3) : Optimize surround cone pooling models using the input cone mosaic STF responses. (NOTE: This step takes a long time)';
    operationDescriptors{6}  = '[ 6] mRGCMosaic visualization       : Inspect optimized cone pooling kernels';
    operationDescriptors{7}  = '[ 7] mRGCMosaic generation (step 4) : Generate the compute-ready mRGCMosaic based on the optimized cone pooling kernels';
    operationDescriptors{8}  = '[ 8] Compute-ready mRGCMosaic validation (step 1) : Compute visual STFs for all cells in the mosaic';
    operationDescriptors{9}  = '[ 9] Compute-ready mRGCMosaic validation (step 2) : Fit the DoG model to the computed visual STFs for all cells in the mosaic';
    operationDescriptors{10} = '[10] Compute-ready mRGCMosaic validation (step 3) : Visualize cone pooling RF maps and visual STF for individual target RGCs';
    operationDescriptors{11} = '[11] Compute-ready mRGCMosaic validation (step 4) : Visualize fitted DoG model params for all cells in the current mosaic';
    operationDescriptors{12} = '[12] Compute-ready mRGCMosaic validation (step 5) : Visualize fitted DoG model params for all cells in multiple mosaics';
    operationDescriptors{13} = '[13] Compute-ready mRGCMosaic (RF computation) : Compute visual RFs (subspace) for all cells in the mosaic';
    operationDescriptors{14} = '[14] Compute-ready mRGCMosaic (RF computation) : Visualize visual RF maps (subspace) for individual target RGCs';

    % 01. Generate the center-connected mosaic
    operationSetToPerformContains.generateCenterConnectedRGCMosaic = ~true;

    % 02. Visualize the center-connected mosaic
    operationSetToPerformContains.visualizeCenterConnectedRGCMosaicAndRemoveUnwantedRGCs = ~true;

    % 03. Visualize PSFs at various positions within the mosaic
    operationSetToPerformContains.visualizePSFsWithinRGCMosaic = ~true;

    % 04 Compute cone mosaic STF responses
    operationSetToPerformContains.computeInputConeMosaicSTFresponses = ~true;

    % 05 Derive optimized surround cone pooling kernels (takes a long time)
    operationSetToPerformContains.optimizeSurroundConePoolingModels = ~true;

    % 06 Inspect optimized cone pooling models
    operationSetToPerformContains.inspectOptimizedSurroundConePoolingModels = ~true;

    % 07 Generate the compute-ready mRGC mosaic
    operationSetToPerformContains.generateComputeReadyMidgetRGCMosaic = ~true;

    % 08 Validation: compute visual STFs for all cells in the mosaic
    operationSetToPerformContains.computeVisualSTFsOfTheComputeReadyMidgetRGCMosaic = ~true;
    
    % 09 Validation: fit the DoG model to the computed visual STFs
    operationSetToPerformContains.fitVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic = ~true;

    % 10 Visualization: Cone pooling RF maps and visual STF 
    operationSetToPerformContains.visualizeConePoolingRFmapAndVisualSTFforTargetRGC = ~true;

    % 11 Visualization: Fitted DoG model params for all cells in the current mosaic
    operationSetToPerformContains.visualizeFittedSTFsAcrossTheComputeReadyMidgetRGCMosaic = ~true;

    % 12 Visualization: Fitted DoG model params for all cells in the current mosaic
    operationSetToPerformContains.visualizeFittedSTFsAcrossMultipleComputeReadyMidgetRGCMosaics = ~true;

    % 13 Subspace RF mapping: compute visual RFs (subspace) for all cells 
    operationSetToPerformContains.computeVisualRFsAcrossTheComputeReadyMidgetRGCMosaic = ~true;

    % 14 Subspace RF mapping: visualize visual RF maps 
    operationSetToPerformContains.visualizeVisualRFmapForTargetRGC = ~true;

    operationSetToPerformContains.animateModelConvergence = ~true;

    invalidActionSelected = true;
    validChoiceIDs = 1:numel(operationDescriptors);

    while (invalidActionSelected)
        % Present options
        fprintf('\n\nAvailable actions for the mosaic at eccentricity (%2.1f,%2.1f):', ...
            mosaicParams.eccDegs(1), mosaicParams.eccDegs(2));
        for iString = 1:numel(operationDescriptors)
            fprintf('\n\t%s', operationDescriptors{iString});
            if (iString == 7) || (iString == 12)
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
                            error('no such field');
                        end

                    case 2
                        % Visualize center-connected mosaic
                        if (isfield(operationSetToPerformContains, 'visualizeCenterConnectedRGCMosaicAndRemoveUnwantedRGCs'))
                            operationSetToPerformContains.visualizeCenterConnectedRGCMosaicAndRemoveUnwantedRGCs = true;
                        else
                            error('no such field');
                        end

                    case 3
                        % Visualize PSFs at various positions within the mosaic
                        if (isfield(operationSetToPerformContains, 'visualizePSFsWithinRGCMosaic'))
                            operationSetToPerformContains.visualizePSFsWithinRGCMosaic = true;
                        else
                            error('no such field');
                        end

                    case 4
                        % Compute input cone mosaic STF responses
                        if (isfield(operationSetToPerformContains, 'computeInputConeMosaicSTFresponses'))
                            operationSetToPerformContains.computeInputConeMosaicSTFresponses = true;
                        else
                            error('no such field');
                        end

                    case 5
                        % Derive optimized surround cone pooling kernels (takes a long time)
                        if (isfield(operationSetToPerformContains, 'optimizeSurroundConePoolingModels'))
                            operationSetToPerformContains.optimizeSurroundConePoolingModels = true;
                        else
                            error('no such field');
                        end

                    case 6
                        % Examine optimized cone pooling models at all or certain  grid nodes
                        if (isfield(operationSetToPerformContains, 'inspectOptimizedSurroundConePoolingModels'))
                            operationSetToPerformContains.inspectOptimizedSurroundConePoolingModels = true;
                        else
                            error('no such field');
                        end

                    case 7
                        % Generate a compute-ready MRGC mosaic
                        if (isfield(operationSetToPerformContains, 'generateComputeReadyMidgetRGCMosaic'))
                            operationSetToPerformContains.generateComputeReadyMidgetRGCMosaic = true;
                        else
                            error('no such field');
                        end

                    case 8
                        % Validate a compute-ready mRGCMosaic: step1 - compute visual STFs for all cells
                        if (isfield(operationSetToPerformContains, 'computeVisualSTFsOfTheComputeReadyMidgetRGCMosaic '))
                            operationSetToPerformContains.computeVisualSTFsOfTheComputeReadyMidgetRGCMosaic  = true;
                        else
                            error('no such field');
                        end
              
                    case 9
                        % Validate a compute-ready mRGCMosaic: step2 - fit a DoG model to all the computed visual STFs for all cells
                        if (isfield(operationSetToPerformContains, 'fitVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic'))
                            operationSetToPerformContains.fitVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic = true;
                        else
                            error('no such field');
                        end

                    case 10
                        % Visualize the derived spatial RF and resulting visual STF for a target RGC 
                        if (isfield(operationSetToPerformContains, 'visualizeConePoolingRFmapAndVisualSTFforTargetRGC'))
                            operationSetToPerformContains.visualizeConePoolingRFmapAndVisualSTFforTargetRGC = true;
                        else
                            error('no such field');
                        end

                    case 11
                        % Validate a compute-ready mRGCMosaic: step3 - visualize fitted visual STFs for all cells in a single mRGC mosaic
                        if (isfield(operationSetToPerformContains, 'visualizeFittedSTFsAcrossTheComputeReadyMidgetRGCMosaic'))
                            operationSetToPerformContains.visualizeFittedSTFsAcrossTheComputeReadyMidgetRGCMosaic = true;
                        else
                            error('no such field');
                        end

                    case 12
                        % Validate a compute-ready mRGCMosaic: step4 - visualize fitted visual STFs for all cells in multiple mRGC mosaics
                        if (isfield(operationSetToPerformContains, 'visualizeFittedSTFsAcrossMultipleComputeReadyMidgetRGCMosaics'))
                            operationSetToPerformContains.visualizeFittedSTFsAcrossMultipleComputeReadyMidgetRGCMosaics = true;
                        else
                            error('no such field');
                        end
                        
                    case 13
                        if (isfield(operationSetToPerformContains, 'computeVisualRFsAcrossTheComputeReadyMidgetRGCMosaic'))
                            operationSetToPerformContains.computeVisualRFsAcrossTheComputeReadyMidgetRGCMosaic = true;
                        else
                            error('no such field');
                        end
                        
                    case 14
                        if (isfield(operationSetToPerformContains, 'visualizeVisualRFmapForTargetRGC '))
                            operationSetToPerformContains.visualizeVisualRFmapForTargetRGC  = true;
                        else
                            error('no such field');
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
