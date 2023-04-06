function performVisualizeConePoolingRFmapAndVisualSTFforTargetRGC(mosaicParams, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('tickSeparationArcMin', 6, @isscalar);
    p.addParameter('visualizedSpatialFrequencyRange', [], @(x)(isempty(x)||(numel(x)==2)));
    p.parse(varargin{:});
    tickSeparationArcMin = p.Results.tickSeparationArcMin;
    visualizedSpatialFrequencyRange = p.Results.visualizedSpatialFrequencyRange;

    % Ask the user which optics were used for computing the input cone
    % mosaic STF responses, so we can obtain the corresponding coneMosaicSTFresponsesFileName
    opticsParams = MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);

    % Ask the user which H1 cell index to use for optimizing the RF
    % surround pooling model
    retinalRFmodelParams = MosaicPoolingOptimizer.chooseRFmodelForSurroundConePoolingOptimization(mosaicParams, opticsParams);

    % Generate the filename of the compute-ready mRGCMosaic to generate
    computeReadyMosaicFileName = MosaicPoolingOptimizer.resourceFileNameAndPath('computeReadyMosaic', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParams, ...
                'retinalRFmodelParams', retinalRFmodelParams);

    % Generate filename for the computed mRGCMosaicSTF responses
    [mRGCMosaicSTFresponsesFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mRGCMosaicSTFresponses', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParams);

    % Get PDF directory
    [~,~,pdfDirectory] = MosaicPoolingOptimizer.resourceFileNameAndPath('pdfsDirectory');

    % Ask the user which RGC to look for:
    % position, # of center cones, majority cone type
    targetRGCposition = input('Enter (xy) position of target RGC (e.g., [5.6 -1.3]): ');
    targetCenterConesNum = input('Enter # of center cones num (e.g, 3): ');
    targetCenterConeMajorityType = input('Enter type of majority center cone num (either cMosaic.LCONE_ID or cMosaic.MCONE_ID): ');

    MosaicPoolingOptimizer.visualizeConePoolingRFmapAndVisualSTFforTargetRGC(...
            fullfile(resourcesDirectory, computeReadyMosaicFileName), ...
            fullfile(resourcesDirectory, mRGCMosaicSTFresponsesFileName), ...
            fullfile(pdfDirectory, 'retinalConePoolingRFmapAndVisualSTF.pdf'), ...
            targetRGCposition, targetCenterConesNum, targetCenterConeMajorityType, ...
            'tickSeparationArcMin', tickSeparationArcMin, ...
            'visualizedSpatialFrequencyRange',  visualizedSpatialFrequencyRange);

end
