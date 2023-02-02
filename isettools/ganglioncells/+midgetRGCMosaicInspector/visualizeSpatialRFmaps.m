function visualizeSpatialRFmaps(mosaicCenterParams, maxRGCsNum)
    
     midgetRGCMosaicInspector.say('Visualizing spatial receptive field maps');

    % Generate mosaic filename and directory
    [mosaicFileName, mosaicDirectory] = midgetRGCMosaicInspector.mosaicFileName(...
        mosaicCenterParams);

    % Load the mosaic
    load(mosaicFileName, 'theMidgetRGCmosaic');

    fprintf(['Center/Surround RFs of theMidgetRGCmosaic were generated based on %d spatial RTVF objects \n' ...
             'for %d different cases of center cones num\n'], ...
             numel(theMidgetRGCmosaic.theRetinaToVisualFieldTransformerOBJList), ...
             numel(unique(theMidgetRGCmosaic.theConesNumPooledByTheRFcenterGrid)));


    % Examine RFs along some meridian angle
    theMeridianAngle = input('Enter meridian angle to visualize : ','s');
    if (isempty(theMeridianAngle))
        theMeridianAngle = 0;
    else
        theMeridianAngle = str2double(theMeridianAngle);
    end

    theMeridianRadius = 1.5 * sqrt(2.0);
    rgcIndicesToAnalyze = midgetRGCMosaicInspector.rgcIndicesAlongMeridianWithAngle(...
            theMeridianAngle, theMeridianRadius, ...
            theMidgetRGCmosaic.rgcRFpositionsDegs, maxRGCsNum);

    videoFileName = sprintf('RFsAlongMeridianAt_%2.1fdegs', theMeridianAngle);
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

        % Retrieve the net weights of the L-cones and the M-cones in the RF center
        [theCenterConeTypeWeights, theCenterConeTypeNum] = theMidgetRGCmosaic.centerConeTypeWeights(theTargetRGCindex);


        % Extract the PSF data to visualize
        if (strcmp(visualizedPSF,'nearest'))
            % Visualize the PSF of the nearest RTVF object (i.e. max weight for this RGC)
            objIndexForVisualizingPSF = nearestRTVFobjectIndex;
        end
        theRTVFTobj = theMidgetRGCmosaic.theRetinaToVisualFieldTransformerOBJList{objIndexForVisualizingPSF};

        theVisualizedPSFData = theRTVFTobj.theSpectrallyWeightedPSFData;
        theVisualizedPSFData.psfSupportXdegs = theVisualizedPSFData.psfSupportXdegs + ...
            theMidgetRGCmosaic.rgcRFpositionsDegs(theTargetRGCindex,1);
        theVisualizedPSFData.psfSupportYdegs = theVisualizedPSFData.psfSupportYdegs + ...
            theMidgetRGCmosaic.rgcRFpositionsDegs(theTargetRGCindex,2);
       

        % Extract the STF of the nearest RTVFobject
        theNearestRTVFTobj = theMidgetRGCmosaic.theRetinaToVisualFieldTransformerOBJList{nearestRTVFobjectIndex};


        fittedVisualSTF = theNearestRTVFTobj.LconeRFcomputeStruct.theSTF.fitted(:);
        theLconeCenterRGCModelVisualSTF = fittedVisualSTF / max(fittedVisualSTF(:));

        fittedVisualSTF = theNearestRTVFTobj.LconeRFcomputeStruct.theSTF.fitted(:);
        theMconeCenterRGCModelVisualSTF = fittedVisualSTF / max(fittedVisualSTF(:));


        spatialFrequencySupport = theNearestRTVFTobj.LconeRFcomputeStruct.theSTF.support(:);
        LconeCenterVisualSTF = theLconeCenterRGCModelVisualSTF * theCenterConeTypeWeights(cMosaic.LCONE_ID)/sum(theCenterConeTypeWeights);
        MconeCenterVisualSTF = theMconeCenterRGCModelVisualSTF * theCenterConeTypeWeights(cMosaic.MCONE_ID)/sum(theCenterConeTypeWeights);
        multiSpectralVisualSTF = theNearestRTVFobjSTFdata.LconeCenterVisualSTF + theNearestRTVFobjSTFdata.MconeCenterVisualSTF;

        [hFig, allAxes] = theMidgetRGCmosaic.visualizeSpatialRFs(...
                'onlyForRGCwithIndex', theTargetRGCindex, ...
                'visualizedRFspatialExtent', 0.3, ...
                'withPSFData', theVisualizedPSFData, ...
                'generateVideo', false, ...
                'withEccentricityCrossHairs', false, ...
                'fontSize', 16);

        % replace graphic is (1,1) with the STFs of the nearest RTVF
        cla(allAxes{1,1});
        plot(allAxes{1,1}, spatialFrequencySupport, LconeCenterVisualSTF , 'r-', 'LineWidth', 1.5);
        hold(allAxes{1,1}, 'on');
        plot(allAxes{1,1}, spatialFrequencySupport, MconeCenterVisualSTF , 'g-', 'LineWidth', 1.5);
        plot(allAxes{1,1}, spatialFrequencySupport, multiSpectralVisualSTF, 'k--', 'LineWidth', 1.5);

        legend(allAxes{1,1},{'nearest L-cone RTVF visual STF', 'nearest M-cone RTVF visual STF', 'nearest mixed cone RTVF visual STF'}, 'Location', 'SouthWest');
        set(allAxes{1,1}, 'YLim', [0 1], 'YTick', 0:0.2:1, 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100], 'XLim', [0.01 100]);
        grid(allAxes{1,1}, 'on')
        xlabel(allAxes{1,1}, 'spatial frequency (c/deg)')
        set(allAxes{1,1}, 'FontSize', 16);

        drawnow;
        videoOBJ.writeVideo(getframe(hFig));
        pause
    end % iRGC

    videoOBJ.close();
    fprintf('spatial RFs video saved at %s\n', videoFileName);
end