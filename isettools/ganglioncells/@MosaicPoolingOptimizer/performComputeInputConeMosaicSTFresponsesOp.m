function performComputeInputConeMosaicSTFresponsesOp(mosaicParams)

    % Generate the mosaic filename
    [mosaicFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mosaic', ...
            'mosaicParams', mosaicParams);

    % Load the generated center-only connected mRGCmosaic
    load(fullfile(resourcesDirectory,mosaicFileName), 'theMidgetRGCMosaic');

    % Ask the user about what optics to use for computing the input cone
    % mosaic STF responses
    [opticsParams, opticsToEmploy, coneMosaicSTFresponsesFileName] = ...
        MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);

    % Generate and set the optics
    theMidgetRGCMosaic.setTheOptics(opticsParams);

    % Instantiate a MosaicPoolingOptimizer object with the center-connected
    % mRGC mosaic and no sampling grid
    theMosaicPoolingOptimizerOBJ = MosaicPoolingOptimizer(...
         theMidgetRGCMosaic, ...
         'generateSamplingGrids', false);


    % Compute the input cone mosaic visual STFs using  
    % using full field stimuli
    stimSizeDegs = [];

    % Positioned at the mosaic's center
    gridNodeIndex = [];

    % Generate and save the input cone mosaic STF responses
    theMosaicPoolingOptimizerOBJ.generateInputConeMosaicSTFresponses(...
            gridNodeIndex, stimSizeDegs, ...
            fullfile(resourcesDirectory, coneMosaicSTFresponsesFileName), ...
            'useParfor', ~true, ...
            'visualizedResponses', ~true, ...
             'opticsToEmploy', opticsToEmploy);


    % Alternatively, we could do this at different grid locations
    compuseConeMosaicSTFResponsesSeparatelyAtEachGridLocation = false;
    if (compuseConeMosaicSTFResponsesSeparatelyAtEachGridLocation)
            % Instantiate the mosaic pooling optimizer
            theMosaicPoolingOptimizerOBJ = MosaicPoolingOptimizer(...
                theMidgetRGCMosaic, ...
                'generateSamplingGrids', true);

            % Compute the input cone mosaic visual STFs, at a single
            % gridNode (using small stimulus patches centered on that node)
            for gridNodeIndex = 1:theMosaicPoolingOptimizerOBJ.gridNodesNum
                % Responses filename
                coneMosaicSTFresponsesFileNameForThisNode = strrep(coneMosaicSTFresponsesFileName, '.mat', sprintf('_AtGridNode_%d.mat', gridNodeIndex));
    
                stimSizeDegs = [1 1];
                theMosaicPoolingOptimizerOBJ.generateInputConeMosaicSTFresponses(...
                    gridNodeIndex, stimSizeDegs,  ...
                    fullfile(resourcesDirectory, coneMosaicSTFresponsesFileNameForThisNode), ...
                    'useParfor', true, ...
                    'visualizedResponses', ~true, ...
                    'opticsToEmploy', opticsToEmploy);
            end
    end

end