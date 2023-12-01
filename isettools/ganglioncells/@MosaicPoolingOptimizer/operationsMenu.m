function operationSetToPerformContains = operationsMenu(mosaicParams)

    operationDescriptors{1}  = '[ 1] mRGCMosaic generation (step G1)     : Generate the center-connected mRGCMosaic ';
    operationDescriptors{2}  = '[ 2] mRGCMosaic visualization            : Visualize the center-connected mRGCMosaic (and possibly remove unwanted RGCs)';
    operationDescriptors{3}  = '[ 3] mRGCMosaic visualization            : Visualize PSFs at a 3x3 grid within the mRGCMosaic';
    operationDescriptors{4}  = '[ 4] mRGCMosaic visualization            : Visualize conePSFs at some location within the mRGCMosaic';

    operationDescriptors{5}  = '[ 5] mRGCMosaic generation (step G2)     : Compute STF responses of the input cone mosaic ';
    operationDescriptors{6}  = '[ 6] inputConeMosaic visualization       : Visualize visually-projected cone Rc';

    operationDescriptors{7}  = '[ 7] mRGCMosaic generation (step G3)     : Optimize surround cone pooling models using the input cone mosaic STF responses. (NOTE: This step takes a long time) ';
    operationDescriptors{8}  = '[ 8] mRGCMosaic visualization            : Inspect optimized cone pooling kernels for a single RGC mosaic';
    operationDescriptors{9}  = '[ 9] mRGCMosaic visualization            : Summarize optimized cone pooling kernels across multiple RGC mosaics';
    operationDescriptors{10}  = '[10] mRGCMosaic generation (step G4)     : Generate the compute-ready mRGCMosaic based on the optimized cone pooling kernels ';
    
    operationDescriptors{11} = '[11] Compute-ready mRGCMosaic validation : Compute visual STFs for all cells in the RGC mosaic';
    operationDescriptors{12} = '[12] Compute-ready mRGCMosaic validation : Fit the DoG model to the computed visual STFs for all cells in the RGC mosaic';
    operationDescriptors{13} = '[13] mRGCMosaic generation (step G5)     : Adjust gain of the compute-ready mRGCMosaic based on the  fitted visual STFs';
    operationDescriptors{14} = '[14] mRGCMosaic generation (step G5)     : Reset gain of the compute-ready mRGCMosaic (1/integrated center weights)';

    operationDescriptors{15} = '[15] Compute-ready mRGCMosaic validation : Visualize cone pooling RF maps and visual STF for individual target RGCs';
    operationDescriptors{16} = '[16] Compute-ready mRGCMosaic validation : Visualize fitted DoG model params for all cells in the RGC mosaic';
    operationDescriptors{17} = '[17] Compute-ready mRGCMosaic validation : Visualize fitted DoG model params for all cells in multiple mRGC mosaics';

    operationDescriptors{18} = '[18] Compute-ready mRGCMosaic (RF computation) : Compute visual RFs (subspace) for cells within the RGC mosaic';
    operationDescriptors{19} = '[19] Compute-ready mRGCMosaic (RF computation) : Compute visual RFs (m-sequence) for cells within the RGC mosaic';


    operationDescriptors{20} = '[20] ComputeVisualRFcenterMapsViaDirectConvolutionWithPSF';

    operationDescriptors{21} = '[21] Compute-ready mRGCMosaic (comparisons) : Contrast STF responses for different optics/stim chromaticities';
    operationDescriptors{22} = '[22] Compute-ready mRGCMosaic (comparisons) : Contrast MSequence RFs for different optics/stim chromaticities';

    % 01. Generate the center-connected mosaic
    operationSetToPerformContains.generateCenterConnectedRGCMosaic = ~true;

    % 02. Visualize the center-connected mosaic
    operationSetToPerformContains.visualizeCenterConnectedRGCMosaicAndRemoveUnwantedRGCs = ~true;

    % 03. Visualize PSFs at various positions within the mosaic
    operationSetToPerformContains.visualizePSFsWithinRGCMosaic = ~true;

    % 04. Visualize cone PSFs at one position within the mosaic
    operationSetToPerformContains.visualizeConePSFsAtLocationWithinRGCMosaic = ~true;

    % 05 Compute cone mosaic STF responses
    operationSetToPerformContains.computeInputConeMosaicSTFresponses = ~true;

    % 06 ComputeInputConeMosaicVisuallyProjectedCharacteristicRadii
    operationSetToPerformContains.computeInputConeMosaicVisuallyProjectedCharacteristicRadii = ~true;

    % 07 Derive optimized surround cone pooling kernels (takes a long time)
    operationSetToPerformContains.optimizeSurroundConePoolingModels = ~true;

    % 08 Inspect optimized cone pooling models for a single mosaic
    operationSetToPerformContains.inspectOptimizedSurroundConePoolingModels = ~true;

    % 09 Summarize optimized cone pooling models across multiple mosaics
    operationSetToPerformContains.summarizeOptimizedSurroundConePoolingModels = ~true;

    % 10 Generate the compute-ready mRGC mosaic
    operationSetToPerformContains.generateComputeReadyMidgetRGCMosaic = ~true;

    % 11 Validation: compute visual STFs for all cells in the mosaic
    operationSetToPerformContains.computeVisualSTFsOfTheComputeReadyMidgetRGCMosaic = ~true;
    
    % 12 Validation: fit the DoG model to the computed visual STFs
    operationSetToPerformContains.fitVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic = ~true;

    % 13 Adjust gain of the compute-ready mRGC mosaic based on fitted visual STFs
    operationSetToPerformContains.adjustGainOfComputeReadyMidgetRGCMosaicBaseOnFittedVisualSTFs = ~true;

    % 14 Reset gain of the compute-ready mRGC mosaic
    operationSetToPerformContains.resetGainOfComputeReadyMidgetRGCMosaic = ~true;

    % 15 Visualization: Cone pooling RF maps and visual STF 
    operationSetToPerformContains.visualizeConePoolingRFmapAndVisualSTFforTargetRGC = ~true;

    % 16 Visualization: Fitted DoG model params for all cells in the current mosaic
    operationSetToPerformContains.visualizeDoGparamsOfVisualSTFsOfSingleMidgetRGCMosaic = ~true;

    % 17 Visualization: Fitted DoG model params for all cells in the current mosaic
    operationSetToPerformContains.visualizeDoGparamsOfVisualSTFsOfMultipleMidgetRGCMosaics = ~true;

    % 18 Visual RF mapping using subspace
    operationSetToPerformContains.computeVisualRFsUsingSubspaceMapping = ~true;

    % 19 Visual RF mapping using m-sequence
    operationSetToPerformContains.computeVisualRFsUsingMSequenceMapping = ~true;

    % 20 Map the visual RF center of mRGCs
    operationSetToPerformContains.computeVisualRFcenterMapsViaDirectConvolutionWithPSF = ~true;

    % 21 Contrast mRGC STFs obtained with different optics and/or different  stimulus chromaticities
    operationSetToPerformContains.contrastSTFsAcrossDifferentOpticsOrChromaticities = ~true;

    % 22. Contrast mSequence RFs obtained with different optics and/or different  stimulus chromaticities
    operationSetToPerformContains.contrastMSequenceRFsAcrossDifferentOpticsOrChromaticities = ~true;


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

            if (iString == 6) || (iString == 9) || (iString == 13) || (iString == 14) || (iString == 17) || (iString == 19) || (iString == 20)
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
                        if (isfield(operationSetToPerformContains, 'visualizeConePSFsAtLocationWithinRGCMosaic'))
                            operationSetToPerformContains.visualizeConePSFsAtLocationWithinRGCMosaic = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'visualizeConePSFsAtLocationWithinRGCMosaic');
                        end

                    case 5
                        % Compute input cone mosaic STF responses
                        if (isfield(operationSetToPerformContains, 'computeInputConeMosaicSTFresponses'))
                            operationSetToPerformContains.computeInputConeMosaicSTFresponses = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'computeInputConeMosaicSTFresponses');
                        end

                    case 6
                        % Compute input cone mosaic visually projected characteristic radii
                        if (isfield(operationSetToPerformContains, 'computeInputConeMosaicVisuallyProjectedCharacteristicRadii'))
                            operationSetToPerformContains.computeInputConeMosaicVisuallyProjectedCharacteristicRadii = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'computeInputConeMosaicVisuallyProjectedCharacteristicRadii');
                        end

                    case 7
                        % Derive optimized surround cone pooling kernels (takes a long time)
                        if (isfield(operationSetToPerformContains, 'optimizeSurroundConePoolingModels'))
                            operationSetToPerformContains.optimizeSurroundConePoolingModels = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'optimizeSurroundConePoolingModels');
                        end

                    case 8
                        % Inspect optimized cone pooling models at all or certain  grid nodes
                        if (isfield(operationSetToPerformContains, 'inspectOptimizedSurroundConePoolingModels'))
                            operationSetToPerformContains.inspectOptimizedSurroundConePoolingModels = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'inspectOptimizedSurroundConePoolingModels');
                        end

                    case 9
                        % Summarize optimized cone pooling models across multiple mosaics at all grid nodes
                        if (isfield(operationSetToPerformContains, 'summarizeOptimizedSurroundConePoolingModels'))
                            operationSetToPerformContains.summarizeOptimizedSurroundConePoolingModels = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'summarizeOptimizedSurroundConePoolingModels');
                        end

                    case 10
                        % Generate a compute-ready MRGC mosaic
                        if (isfield(operationSetToPerformContains, 'generateComputeReadyMidgetRGCMosaic'))
                            operationSetToPerformContains.generateComputeReadyMidgetRGCMosaic = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'generateComputeReadyMidgetRGCMosaic');
                        end

                    case 11
                        % Validate a compute-ready mRGCMosaic: step1 - compute visual STFs for all cells
                        if (isfield(operationSetToPerformContains, 'computeVisualSTFsOfTheComputeReadyMidgetRGCMosaic'))
                            operationSetToPerformContains.computeVisualSTFsOfTheComputeReadyMidgetRGCMosaic  = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'computeVisualSTFsOfTheComputeReadyMidgetRGCMosaic');
                        end
              
                    case 12
                        % Validate a compute-ready mRGCMosaic: step2 - fit a DoG model to all the computed visual STFs for all cells
                        if (isfield(operationSetToPerformContains, 'fitVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic'))
                            operationSetToPerformContains.fitVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'fitVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic');
                        end

                    case 13
                        % Adjust gain of compute-read MRGCMosaic based on the computed visual STFs 
                        if (isfield(operationSetToPerformContains, 'adjustGainOfComputeReadyMidgetRGCMosaicBaseOnFittedVisualSTFs'))
                            operationSetToPerformContains.adjustGainOfComputeReadyMidgetRGCMosaicBaseOnFittedVisualSTFs = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'adjustGainOfComputeReadyMidgetRGCMosaicBaseOnFittedVisualSTFs');
                        end

                    case 14
                        % Reset gain of compute-read MRGCMosaic 
                        if (isfield(operationSetToPerformContains, 'resetGainOfComputeReadyMidgetRGCMosaic'))
                            operationSetToPerformContains.resetGainOfComputeReadyMidgetRGCMosaic = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'resetGainOfComputeReadyMidgetRGCMosaic');
                        end

                    case 15
                        % Visualize the derived spatial RF and resulting visual STF for a target RGC 
                        if (isfield(operationSetToPerformContains, 'visualizeConePoolingRFmapAndVisualSTFforTargetRGC'))
                            operationSetToPerformContains.visualizeConePoolingRFmapAndVisualSTFforTargetRGC = true;
                        else
                             error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'visualizeConePoolingRFmapAndVisualSTFforTargetRGC');
                        end

                    case 16
                        % Validate a compute-ready mRGCMosaic: step3 - visualize fitted visual STFs for all cells in a single mRGC mosaic
                        if (isfield(operationSetToPerformContains, 'visualizeDoGparamsOfVisualSTFsOfSingleMidgetRGCMosaic'))
                            operationSetToPerformContains.visualizeDoGparamsOfVisualSTFsOfSingleMidgetRGCMosaic = true;
                        else
                           error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'visualizeDoGparamsOfVisualSTFsOfSingleMidgetRGCMosaic');
                        end

                    case 17
                        % Validate a compute-ready mRGCMosaic: step4 - visualize fitted visual STFs for all cells in multiple mRGC mosaics
                        if (isfield(operationSetToPerformContains, 'visualizeDoGparamsOfVisualSTFsOfMultipleMidgetRGCMosaics'))
                            operationSetToPerformContains.visualizeDoGparamsOfVisualSTFsOfMultipleMidgetRGCMosaics = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'visualizeDoGparamsOfVisualSTFsOfMultipleMidgetRGCMosaics');
                        end
                        
                    case 18
                        if (isfield(operationSetToPerformContains, 'computeVisualRFsUsingSubspaceMapping'))
                            operationSetToPerformContains.computeVisualRFsUsingSubspaceMapping = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: %s''.', 'computeVisualRFsUsingSubspaceMapping');
                        end
                        
                    case 19
                        if (isfield(operationSetToPerformContains, 'computeVisualRFsUsingMSequenceMapping'))
                            operationSetToPerformContains.computeVisualRFsUsingMSequenceMapping  = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: ''%s''.', 'computeVisualRFsUsingMSequenceMapping');
                        end
                       
                    case 20
                        if (isfield(operationSetToPerformContains, 'computeVisualRFcenterMapsViaDirectConvolutionWithPSF'))
                            operationSetToPerformContains.computeVisualRFcenterMapsViaDirectConvolutionWithPSF = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: ''%s''.', 'computeVisualRFcenterMapsViaDirectConvolutionWithPSF');
                        end

                    case 21
                        if (isfield(operationSetToPerformContains, 'contrastSTFsAcrossDifferentOpticsOrChromaticities'))
                            operationSetToPerformContains.contrastSTFsAcrossDifferentOpticsOrChromaticities = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: ''%s''.', 'contrastSTFsAcrossDifferentOpticsOrChromaticities');
                        end

                    case 22
                        if (isfield(operationSetToPerformContains, 'contrastMSequenceRFsAcrossDifferentOpticsOrChromaticities'))
                            operationSetToPerformContains.contrastMSequenceRFsAcrossDifferentOpticsOrChromaticities = true;
                        else
                            error('MosaicPoolingOptimizer.operationsMenu: no such field: ''%s''.', 'contrastMSequenceRFsAcrossDifferentOpticsOrChromaticities');
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
