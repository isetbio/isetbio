function performVisualizeConePoolingRFmapAndVisualSTFforTargetRGC(mosaicParams, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('tickSeparationArcMin', 6, @isscalar);
    p.addParameter('normalizedPeakSurroundSensitivity', 0.4, @isscalar);
    p.addParameter('visualizedSpatialFrequencyRange', [], @(x)(isempty(x)||(numel(x)==2)));
    p.addParameter('reverseXDir', false, @islogical);
    p.addParameter('gridlessLineWeightingFunctions', false, @islogical);

    p.parse(varargin{:});
    tickSeparationArcMin = p.Results.tickSeparationArcMin;
    visualizedSpatialFrequencyRange = p.Results.visualizedSpatialFrequencyRange;
    normalizedPeakSurroundSensitivity = p.Results.normalizedPeakSurroundSensitivity;
    reverseXDir = p.Results.reverseXDir;
    gridlessLineWeightingFunctions = p.Results.gridlessLineWeightingFunctions;

    % Ask the user what optics were used for computing the compute-ready MRGC mosaic
    fprintf('\n---> Select the optics that were used to compute the compute-ready mosaic\n');
    opticsParamsForComputeReadyMosaic = ...
        MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);

    % Ask the user which H1 cell index to use for optimizing the RF
    % surround pooling model
    retinalRFmodelParams = MosaicPoolingOptimizer.chooseRFmodelForSurroundConePoolingOptimization(...
        mosaicParams, opticsParamsForComputeReadyMosaic);

    % Generate the filename of the compute-ready mRGCMosaic to use
    [computeReadyMosaicFileName, computeReadyMosaicResourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('computeReadyMosaic', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParamsForComputeReadyMosaic, ...
                'retinalRFmodelParams', retinalRFmodelParams);


    % Now, ask the user what optics were used for computing the input cone
    % mosaic STF responses, so we can obtain the corresponding coneMosaicSTFresponsesFileName
    fprintf('\n---> Select the optics that were used to compute the input cone mosaic STF responses on which the mRGC mosaic STF responses were based on\n');
    opticsParamsForMRGCSTFs = ...
        MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);

    % Generate filename for the computed mRGCMosaicSTF responses
    [mRGCMosaicSTFresponsesFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mRGCMosaicSTFresponses', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParamsForMRGCSTFs);

    fprintf('\n---> Select the chromaticity that was used to compute the input cone mosaic STF responses on which the mRGC mosaic STF responses were based on\n');
    % Ask the user what stimulus chromaticity to use
    [~, ~, mRGCMosaicSTFresponsesFileName] = ...
        MosaicPoolingOptimizer.chooseStimulusChromaticityForMosaicResponsesAndUpdateFileName(...
            mRGCMosaicSTFresponsesFileName, 'STFresponses');

    
    % Get PDF directory
    [~,~,pdfDirectory] = MosaicPoolingOptimizer.resourceFileNameAndPath('pdfsDirectory', ...
        'mosaicParams', mosaicParams);

    % Ask the user which RGC to look for:
    % position, # of center cones, majority cone type
    rgcSpecification = input('Plot RGC with specific index (1), or RGC at a target position (2) ? ');
    if (rgcSpecification == 1)
        targetRGCposition = [];
        targetCenterConesNum = [];
        targetCenterConeMajorityType = [];
    else
        targetRGCposition = input('Enter (xy) position of target RGC (e.g., [5.6 -1.3]): ');
        targetCenterConesNum = input('Enter # of center cones num (e.g, 3): ');
        targetCenterConeMajorityType = input('Enter type of majority center cone num (either cMosaic.LCONE_ID or cMosaic.MCONE_ID): ');
    end

    MosaicPoolingOptimizer.visualizeConePoolingRFmapAndVisualSTFforTargetRGC(...
            fullfile(computeReadyMosaicResourcesDirectory, computeReadyMosaicFileName), ...
            fullfile(resourcesDirectory, mRGCMosaicSTFresponsesFileName), ...
            fullfile(pdfDirectory, 'retinalConePoolingRFmapAndVisualSTF.pdf'), ...
            targetRGCposition, targetCenterConesNum, targetCenterConeMajorityType, ...
            'tickSeparationArcMin', tickSeparationArcMin, ...
            'normalizedPeakSurroundSensitivity', normalizedPeakSurroundSensitivity, ...
            'visualizedSpatialFrequencyRange',  visualizedSpatialFrequencyRange, ...
            'reverseXDir', reverseXDir, ...
            'gridlessLineWeightingFunctions', gridlessLineWeightingFunctions);

end
