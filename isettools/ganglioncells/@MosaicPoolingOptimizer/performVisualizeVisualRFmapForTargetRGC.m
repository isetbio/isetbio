function performVisualizeVisualRFmapForTargetRGC(mosaicParams, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('tickSeparationArcMin', 6, @isscalar);
    p.addParameter('reverseXDir', false, @islogical);
    p.addParameter('gridlessLineWeightingFuncions', false, @islogical);
    p.parse(varargin{:});
    tickSeparationArcMin = p.Results.tickSeparationArcMin;
    reverseXDir = p.Results.reverseXDir;
    gridlessLineWeightingFuncions = p.Results.gridlessLineWeightingFuncions;

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
    [~, mRGCMosaicSubspaceResponsesFileName] = ...
        MosaicPoolingOptimizer.chooseStimulusChromaticityForMosaicResponsesAndUpdateFileName(...
        mRGCMosaicSubspaceResponsesFileName, 'SubspaceResponses');

    % Optimally generated RF maps filename
    optimallyMappedSubspaceRFmapsFileName = strrep(mRGCMosaicSubspaceResponsesFileName, '.mat', '_optimallyMappedRFs.mat');

    % PDF directory
    [~,~,pdfDirectory] = MosaicPoolingOptimizer.resourceFileNameAndPath('pdfsDirectory', ...
        'mosaicParams', mosaicParams);

    MosaicPoolingOptimizer.visualizeVisualRFmapsForMultipleTargetRGCs(...
            theComputeReadyMRGCmosaic, ...
            fullfile(resourcesDirectory, optimallyMappedSubspaceRFmapsFileName), ...
            fullfile(resourcesDirectory, mRGCMosaicSubspaceResponsesFileName), ...
            fullfile(pdfDirectory, 'visualRFmap.pdf'), ...
            'tickSeparationArcMin', tickSeparationArcMin, ...
            'reverseXDir', reverseXDir, ...
            'gridlessLineWeightingFuncions', gridlessLineWeightingFuncions);
end
