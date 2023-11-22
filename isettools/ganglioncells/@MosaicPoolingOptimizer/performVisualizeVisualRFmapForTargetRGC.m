function performVisualizeVisualRFmapForTargetRGC(mosaicParams, stimPositionDegs, rfMappingPixelMagnificationFactor, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('tickSeparationArcMin', 6, @isscalar);
    p.addParameter('reverseXDir', false, @islogical);
    p.addParameter('gridlessLineWeightingFunctions', false, @islogical);
    p.addParameter('maxSFLimit', [], @(x)(isempty(x)||isscalar(x)));
    p.addParameter('maxSFToBeAnalyzed', [], @(x)(isempty(x)||isscalar(x)));
    p.parse(varargin{:});
    tickSeparationArcMin = p.Results.tickSeparationArcMin;
    reverseXDir = p.Results.reverseXDir;
    gridlessLineWeightingFunctions = p.Results.gridlessLineWeightingFunctions;

     % stimulus max spatial frequency (can be set to lower that optimal retinal sf)
    maxSFLimit = p.Results.maxSFLimit;

    % max spatial frequency to include in the RFmap
    maxSFToBeAnalyzed = p.Results.maxSFToBeAnalyzed;
    
    % Ask the user which optics were used for computing the input cone
    % mosaic STF responses, so we can obtain the corresponding coneMosaicSTFresponsesFileName
    fprintf('\n---> Select the optics that were used to compute the compute-ready mosaic\n');
    opticsParamsForComputeReadyMosaic = MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);

    % Ask the user which H1 cell index to use for optimizing the RF
    % surround pooling model
    retinalRFmodelParams = MosaicPoolingOptimizer.chooseRFmodelForSurroundConePoolingOptimization(mosaicParams, opticsParamsForComputeReadyMosaic);

    % Generate the filename of the compute-ready mRGCMosaic to generate
    [computeReadyMosaicFileName, computeReadyMosaicResourcesDirectory] = ...
            MosaicPoolingOptimizer.resourceFileNameAndPath('computeReadyMosaic', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParamsForComputeReadyMosaic, ...
                'retinalRFmodelParams', retinalRFmodelParams);

    % Load the compute-ready MRGC mosaic
    load(fullfile(computeReadyMosaicResourcesDirectory, computeReadyMosaicFileName), 'theComputeReadyMRGCmosaic');

    % Now, ask the user what optics were used for computing the input cone
    % mosaic STF responses, so we can obtain the corresponding coneMosaicSTFresponsesFileName
    fprintf('\n---> Select the optics that were used to compute the mRGCMosaicSubspaceResponses\n');
    opticsParamsForMRGCSubSpaceResponses = ...
        MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);

    % Generate filename for the mRGCMosaic subspace responses
    [mRGCMosaicSubspaceResponsesFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mRGCMosaicSubspaceResponses', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParamsForMRGCSubSpaceResponses);

    fprintf('\n---> Select the chromaticity that was used to compute the mRGCMosaicSubspaceResponses\n');
    % Ask the user what stimulus chromaticity to use
    [stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, mRGCMosaicSubspaceResponsesFileName] = ...
        MosaicPoolingOptimizer.chooseStimulusChromaticityForMosaicResponsesAndUpdateFileName(...
        mRGCMosaicSubspaceResponsesFileName, 'mRGCMosaicSubspaceResponses');

    % % Optimally generated RF maps filename
    % optimallyMappedSubspaceRFmapsFileName = strrep(mRGCMosaicSubspaceResponsesFileName, '.mat', '_optimallyMappedRFs.mat');
    % 
    % % Add the position info
    % positionPostFix = sprintf('_atPosition_%2.2f_%2.2f_PixelMagnification_%2.3f.mat', stimPositionDegs(1), stimPositionDegs(2), rfMappingPixelMagnificationFactor);
    % optimallyMappedSubspaceRFmapsFileName = strrep(optimallyMappedSubspaceRFmapsFileName, '.mat', positionPostFix);

    if (isempty(maxSFToBeAnalyzed)) || (maxSFToBeAnalyzed>maxSFLimit)
        maxSFToBeAnalyzed = maxSFLimit;
    end

    % Optimally generated RF maps filename
    optimallyMappedSubspaceRFmapsFileName = strrep(mRGCMosaicSubspaceResponsesFileName, '.mat', '_optimallyMappedRFs.mat');

    if (~isempty(maxSFToBeAnalyzed))
        optimallyMappedSubspaceRFmapsFileName = strrep(optimallyMappedSubspaceRFmapsFileName, '.mat', sprintf('_LimitedTo%2.2fCPD.mat', maxSFToBeAnalyzed));
    end

    % PDF directory
    [~,~,pdfDirectory] = MosaicPoolingOptimizer.resourceFileNameAndPath('pdfsDirectory', ...
        'mosaicParams', mosaicParams);

    MosaicPoolingOptimizer.visualizeVisualRFmapsForMultipleTargetRGCs(...
            theComputeReadyMRGCmosaic, ...
            fullfile(resourcesDirectory, optimallyMappedSubspaceRFmapsFileName), ...
            fullfile(pdfDirectory, 'visualRFmap.pdf'), ...
            'tickSeparationArcMin', tickSeparationArcMin, ...
            'reverseXDir', reverseXDir, ...
            'gridlessLineWeightingFunctions', gridlessLineWeightingFunctions);
end
