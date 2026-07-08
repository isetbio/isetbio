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

    % Set the peak gains based on the model cell's integrated center cone
    % weights. This is the default and is what is used by the
    % mRGCMosaic.compute() method, if the mRGCMosaic.rgcRFgains property is
    % empty (not set).
    method = '1/integrated center cone weights';
    maxGain = 1.0;
    methodParams = maxGain;
    theComputeReadyMRGCmosaic.setPeakGains(method, methodParams);


    useVisualRcEstimates = false;
    if (useVisualRcEstimates)

        % Generate filename for the computed mRGCMosaicSTF responses
        [mRGCMosaicSTFresponsesFileName, resourcesDirectory] = ...
            MosaicPoolingOptimizer.resourceFileNameAndPath('mRGCMosaicSTFresponses', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParams);
        mRGCMosaicSTFresponsesFileName = fullfile(resourcesDirectory, mRGCMosaicSTFresponsesFileName);
        
        % Check that the responses filename contain visualSTF data
        variableInfo = who('-file', mRGCMosaicSTFresponsesFileName);
        if (~ismember('theMRGCMosaicOptimalSTFs', variableInfo))
            error('Did not find fitted visual STFs in the specified responses filename (%s).', ...
                mRGCMosaicSTFresponsesFileName);
        end
    
    
        % Load the fitted optimal STFs which contain the visualRcDegs obtained by
        % fitting the composite visual STF. Also load the visualRcDegs estimate
        % obtained by fitting the center visual STF.
        load(mRGCMosaicSTFresponsesFileName, ...
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
        foveolarVisualRcDegs = 0.456/60; %min(visualRcDegsEstimateFromFittingCompositeSTF)
        
        methodParams(numel(methodParams)+1) = foveolarVisualRcDegs;
        theComputeReadyMRGCmosaic.setPeakGains(method, methodParams);
    end


    debugGains = false;
    if (debugGains)
        % Display the range of rgcRFgains
        rgcRFgainsRange = [min(theComputeReadyMRGCmosaic.rgcRFgains(:)) max(theComputeReadyMRGCmosaic.rgcRFgains(:))]
    
        % Set the peak gains based on the model cell's integrated retinal cone aperture areas
        method = '1/integrated center retinal cone apertures'
        % Gain selected to match the peak gain computed by 'CK95formulaAppliedToMRGCMosaicVisualRcDegs' at the fovea
        maxGain = 1.2024e-12;
        methodParams = maxGain;
        theComputeReadyMRGCmosaic.setPeakGains(method, methodParams);
    
        % Display the range of rgcRFgains
        rgcRFgainsRange = [min(theComputeReadyMRGCmosaic.rgcRFgains(:)) max(theComputeReadyMRGCmosaic.rgcRFgains(:))]
    end

    % Save the updated mosaic
    save(computeReadyMosaicFileName, 'theComputeReadyMRGCmosaic', '-v7.3');
    fprintf('The updated rgcRFgains - compute-ready mRGCMosaic was exported to %s overwriting the previous compute-ready mosaic.\n', computeReadyMosaicFileName);

end
