function visualizeSpatialRFmaps(mosaicCenterParams, maxRGCsNum)
    
    % Generate mosaic filename and directory
    [mosaicFileName, mosaicDirectory] = midgetRGCMosaicInspector.mosaicFileName(...
        mosaicCenterParams);

    % RGC indices to visualize
    rgcIndicesToAnalyze = midgetRGCMosaicInspector.rgcIndicesToAnalyze(...
        mosaicFileName, ...
        'maxRGCsNum', maxRGCsNum);

    % Load the mosaic
    load(mosaicFileName, 'theMidgetRGCmosaic');

    fprintf(['Center/Surround RFs of theMidgetRGCmosaic were generated based on %d RTVF objects\n' ...
             'for %d different cases of center cones num\n'], ...
             numel(theMidgetRGCmosaic.theRetinaToVisualFieldTransformerOBJList), ...
             numel(unique(theMidgetRGCmosaic.theConesNumPooledByTheRFcenterGrid)));

    videoFileName = 'RFs';
    videoOBJ = VideoWriter(videoFileName, 'MPEG-4');
    videoOBJ.FrameRate = 10;
    videoOBJ.Quality = 100;
    videoOBJ.open();

    
    % Visualize the PSF of the RTVF object that was closest to a particular mosaic position
    targetMosaicPosition = [0 0];
    [~, objIndexForVisualizingPSF] = min(sum((bsxfun(@minus, theMidgetRGCmosaic.theSamplingPositionGrid, targetMosaicPosition)).^2,2));
    visualizedPSF = 'target positions';

    % Or visualize the PSF of the RTVF object that is closest to each RGC
    visualizedPSF = 'nearest';


    for iRGC = 1:numel(rgcIndicesToAnalyze)
        % Retrieve theTargetRGCindex
        theTargetRGCindex = rgcIndicesToAnalyze(iRGC);

        [triangulatingRTVFobjIndices, triangulatingRTVFobjWeights] = ...
             theMidgetRGCmosaic.triangulatingRTVFobjectIndicesAndWeights(theTargetRGCindex);
        [~,idx] = max(triangulatingRTVFobjWeights);
        nearestRTVFobjectIndex = triangulatingRTVFobjIndices(idx);

        % Extract the PSF data to visualize
        if (strcmp(visualizedPSF,'nearest'))
            % Visualize the PSF of the nearest RTVF object (i.e. max weight for this RGC)
            objIndexForVisualizingPSF = nearestRTVFobjectIndex;
        end
        theRTVFTobj = theMidgetRGCmosaic.theRetinaToVisualFieldTransformerOBJList{objIndexForVisualizingPSF};
        theVisualizedPSFData = theRTVFTobj.theVlambdaWeightedPSFData;
        theVisualizedPSFData.psfSupportXdegs = theVisualizedPSFData.psfSupportXdegs + ...
            theMidgetRGCmosaic.rgcRFpositionsDegs(theTargetRGCindex,1);
        theVisualizedPSFData.psfSupportYdegs = theVisualizedPSFData.psfSupportYdegs + ...
            theMidgetRGCmosaic.rgcRFpositionsDegs(theTargetRGCindex,2);
       

        % Extract the STF of the nearest RTVFobject
        theNearestRTVFTobj = theMidgetRGCmosaic.theRetinaToVisualFieldTransformerOBJList{nearestRTVFobjectIndex};

        % Normalize theRTVFTobj.fitted to the theRTVFTobj.targetSTF
        scaleFactorBasedOnLowestSF = theNearestRTVFTobj.rfComputeStruct.theSTF.target(1)/theNearestRTVFTobj.rfComputeStruct.theSTF.fitted(1);
 
        % Assemble the nearest RTVF object STF data
        theNearestRTVFobjSTFdata = struct();
        theNearestRTVFobjSTFdata.spatialFrequencySupport = theNearestRTVFTobj.rfComputeStruct.theSTF.support(:);
        % The model-achieved STF
        theNearestRTVFobjSTFdata.fittedSTF = theNearestRTVFTobj.rfComputeStruct.theSTF.fitted(:)*scaleFactorBasedOnLowestSF;
        % The model-target STF
        theNearestRTVFobjSTFdata.targetSTF = theNearestRTVFTobj.rfComputeStruct.theSTF.target(:);


        [hFig, allAxes] = theMidgetRGCmosaic.visualizeSpatialRFs(...
                'onlyForRGCwithIndex', theTargetRGCindex, ...
                'visualizedRFspatialExtent', 0.3, ...
                'withPSFData', theVisualizedPSFData, ...
                'generateVideo', false, ...
                'withEccentricityCrossHairs', true, ...
                'fontSize', 16);

        % replace graphic is (1,1) with the STFs of the nearest RTVF
        cla(allAxes{1,1});
        plot(allAxes{1,1}, theNearestRTVFobjSTFdata.spatialFrequencySupport, theNearestRTVFobjSTFdata.targetSTF, 'k-', 'LineWidth', 1.5);
        hold(allAxes{1,1}, 'on');
        plot(allAxes{1,1}, theNearestRTVFobjSTFdata.spatialFrequencySupport, theNearestRTVFobjSTFdata.fittedSTF, 'r-', 'LineWidth', 1.5);
        legend(allAxes{1,1},{'RTVF: target', 'RTVF: fitted'}, 'Location', 'SouthWest');
        set(allAxes{1,1}, 'YLim', [0 1], 'YTick', 0:0.2:1, 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100], 'XLim', [0.01 100]);
        grid(allAxes{1,1}, 'on')
        xlabel(allAxes{1,1}, 'spatial frequency (c/deg)')
        set(allAxes{1,1}, 'FontSize', 16);

        drawnow;
        videoOBJ.writeVideo(getframe(hFig));
    end % iRGC

    videoOBJ.close();
    fprintf('spatial RFs video saved at %s\n', videoFileName);
end