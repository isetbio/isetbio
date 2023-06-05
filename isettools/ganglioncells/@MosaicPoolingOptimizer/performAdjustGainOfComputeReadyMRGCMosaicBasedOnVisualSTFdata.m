function performAdjustGainOfComputeReadyMRGCMosaicBasedOnVisualSTFdata(mosaicParams)

    % Ask the user which optics were used for computing the input cone
    % mosaic STF responses, so we can obtain the corresponding coneMosaicSTFresponsesFileName
    opticsParams = MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);

    % Ask the user which H1 cell index to use for optimizing the RF
    % surround pooling model
    retinalRFmodelParams = MosaicPoolingOptimizer.chooseRFmodelForSurroundConePoolingOptimization(mosaicParams, opticsParams);

    % Generate the filename of the compute-ready mRGCMosaic to generate
    [computeReadyMosaicFileName, computeReadyMosaicResourcesDirectory] = MosaicPoolingOptimizer.resourceFileNameAndPath('computeReadyMosaic', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParams, ...
                'retinalRFmodelParams', retinalRFmodelParams);

    % Load the compute-ready MRGC mosaic
    computeReadyMosaicFileName = fullfile(computeReadyMosaicResourcesDirectory, computeReadyMosaicFileName);
    load(computeReadyMosaicFileName, 'theComputeReadyMRGCmosaic');

    % Generate filename for the computed mRGCMosaicSTF responses
    [mRGCMosaicSTFresponsesFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mRGCMosaicSTFresponses', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParams);

    % Load the STFresponsesFilename
    mRGCMosaicSTFresponsesFileName = fullfile(resourcesDirectory, mRGCMosaicSTFresponsesFileName);
    load( mRGCMosaicSTFresponsesFileName, ...
        'visualRcDegsEstimates', ...
        'theMRGCMosaicOptimalSTFs');

    % For each cell, extract the visualRcDegs estimates obtained by fitting
    % the centerSTF alone and by fitting the compositeSTF
    rgcsNum = numel(theMRGCMosaicOptimalSTFs);
    visualRcDegsEstimateFromFittingCompositeSTF = zeros(1, rgcsNum);
    
    % For each cell extract the visualRFdegs obtained by fitting the composite STF at the
    % optimal orientation
    for iRGC = 1:theComputeReadyMRGCmosaic.rgcsNum

        % Retrieve the fitted STFdata
        theSTFdata = theMRGCMosaicOptimalSTFs{iRGC};
        if (iRGC == 1)
            % Retrieve the index of the RcDegs param
            idxRcDegs = find(strcmp(theSTFdata.DoGfitParams.names, 'RcDegs'));
        end
    
        % Get it
        visualRcDegsEstimateFromFittingCompositeSTF(iRGC) = theSTFdata.DoGfitParams.finalValues(idxRcDegs);
    end

    % Set the peak gains based on the model cell's computed visual Rc Degs and
    % the C&K '95 formula
    method = 'CK95formulaAppliedToMRGCMosaicVisualRcDegs';
    methodParams = visualRcDegsEstimateFromFittingCompositeSTF;
    theComputeReadyMRGCmosaic.setPeakGains(method, methodParams);

    debugGains = false;
    if (debugGains)
        % Display the range of rgcRFgains
        rgcRFgainsRange = [min(theComputeReadyMRGCmosaic.rgcRFgains(:)) max(theComputeReadyMRGCmosaic.rgcRFgains(:))]
    
    
        % Set the peak gains based on the model cell's integrated center cone weights
        method = '1/integrated center cone weights';
        % Gain selected to match the peak gain computed by 'CK95formulaAppliedToMRGCMosaicVisualRcDegs' at the fovea
        maxGain = 3.3*1e3;
        methodParams = maxGain;
        theComputeReadyMRGCmosaic.setPeakGains(method, methodParams);
    
        % Display the range of rgcRFgains
        rgcRFgainsRange = [min(theComputeReadyMRGCmosaic.rgcRFgains(:)) max(theComputeReadyMRGCmosaic.rgcRFgains(:))]
    
        % Set the peak gains based on the model cell's integrated retinal cone aperture areas
        method = '1/integrated center retinal cone apertures';
        % Gain selected to match the peak gain computed by 'CK95formulaAppliedToMRGCMosaicVisualRcDegs' at the fovea
        maxGain = 4e-9;
        methodParams = maxGain;
        theComputeReadyMRGCmosaic.setPeakGains(method, methodParams);
    
        % Display the range of rgcRFgains
        rgcRFgainsRange = [min(theComputeReadyMRGCmosaic.rgcRFgains(:)) max(theComputeReadyMRGCmosaic.rgcRFgains(:))]
    end

    % Save the updated mosaic
    save(computeReadyMosaicFilename, 'theComputeReadyMRGCmosaic', '-v7.3');
    fprintf('The updated rgcRFgains - compute-ready mRGCMosaic was exported in %s.\n', computeReadyMosaicFilename);

end
