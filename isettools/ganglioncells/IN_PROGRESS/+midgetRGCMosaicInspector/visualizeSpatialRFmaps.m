function visualizeSpatialRFmaps(mosaicCenterParams, maxRGCsNum)
    
     midgetRGCMosaicInspector.say('Visualizing spatial receptive field maps');

    % Generate mosaic filename and directory
    [mosaicFileName, mosaicDirectory] = midgetRGCMosaicInspector.mosaicFileName(...
        mosaicCenterParams);

    % Load the mosaic
    load(mosaicFileName, 'theMidgetRGCmosaic');

    fprintf(['Center/Surround RFs of theMidgetRGCmosaic were generated based on %d spatial RTVF objects \n' ...
             'for %d different cases of center cones num\n'], ...
             numel(theMidgetRGCmosaic.theRTVFobjList), ...
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


    for iRGC = 1:numel(rgcIndicesToAnalyze)
        % Retrieve theTargetRGCindex
        theTargetRGCindex = rgcIndicesToAnalyze(iRGC);

        [triangulatingRTVFobjIndices, triangulatingRTVFobjWeights] = ...
             theMidgetRGCmosaic.triangulatingRTVFobjectIndicesAndWeights(theTargetRGCindex);
        [~,idx] = max(triangulatingRTVFobjWeights);
        nearestRTVFobjectIndex = triangulatingRTVFobjIndices(idx);
        theRemainingNearestRTVFobjectIndices = setdiff(triangulatingRTVFobjIndices, nearestRTVFobjectIndex);

        % Retrieve the net weights of the L-cones and the M-cones in the RF center
        [theCenterConeTypeWeights, theCenterConeTypeNum] = theMidgetRGCmosaic.centerConeTypeWeights(theTargetRGCindex);


        % Extract the PSF data to visualize      
        theRTVFTobj = theMidgetRGCmosaic.theRTVFobjList{nearestRTVFobjectIndex};

        theVisualizedPSFData = theRTVFTobj.spectrallyWeightedPSFData;
        theVisualizedPSFData.psfSupportXdegs = theVisualizedPSFData.psfSupportXdegs + ...
            theMidgetRGCmosaic.rgcRFpositionsDegs(theTargetRGCindex,1);
        theVisualizedPSFData.psfSupportYdegs = theVisualizedPSFData.psfSupportYdegs + ...
            theMidgetRGCmosaic.rgcRFpositionsDegs(theTargetRGCindex,2);
       

        % Extract the STF of the nearest RTVFobject
        theNearestRTVFTobj = theMidgetRGCmosaic.theRTVFobjList{nearestRTVFobjectIndex};

        [spatialFrequencySupport, ...
            LconeCenterVisualSTF, MconeCenterVisualSTF, ...
            multiSpectralVisualSTF] = extractPSFsFromRTVFobj(theNearestRTVFTobj, theCenterConeTypeWeights);

        if (~isempty(theRemainingNearestRTVFobjectIndices))
            RTVFobj2 = theMidgetRGCmosaic.theRTVFobjList{theRemainingNearestRTVFobjectIndices(1)};
            RTVFobj3 = theMidgetRGCmosaic.theRTVFobjList{theRemainingNearestRTVFobjectIndices(2)};
            [~, ~, ~, multiSpectralVisualSTF2] = extractPSFsFromRTVFobj(RTVFobj2, theCenterConeTypeWeights);
            [~, ~, ~, multiSpectralVisualSTF3] = extractPSFsFromRTVFobj(RTVFobj3, theCenterConeTypeWeights);
        end


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
        if (~isempty(theRemainingNearestRTVFobjectIndices))
            plot(allAxes{1,1}, spatialFrequencySupport, multiSpectralVisualSTF2, 'k--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
            plot(allAxes{1,1}, spatialFrequencySupport, multiSpectralVisualSTF3, 'k--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
            
            legend(allAxes{1,1},{...
                'nearest L-cone RTVF visual STF', ...
                'nearest M-cone RTVF visual STF', ...
                'nearest RTVF visual STF', ...
                'RTVF #2 visual STF', 'RTVF #3 visual STF'}, 'Location', 'SouthWest');
        else
            legend(allAxes{1,1},{...
                'nearest L-cone RTVF visual STF', ...
                'nearest M-cone RTVF visual STF', ...
                'nearest RTVF visual STF'}, 'Location', 'SouthWest');
        end

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


function  [spatialFrequencySupport, LconeCenterVisualSTF, MconeCenterVisualSTF, multiSpectralVisualSTF] =  ...
    extractPSFsFromRTVFobj(theRTVFTobj, theCenterConeTypeWeights)

    % Compute the L-cone visual RF map
    theLconeVisualRFmap = theRTVFTobj.visualRFfromRetinalConePooling(...
                theRTVFTobj.LconeRFcomputeStruct.modelConstants, ...
                theRTVFTobj.LconeRFcomputeStruct.retinalConePoolingParams.finalValues);
  
    % Compute the L-cone visual STF data for the L-cone visual RF map
    theLconeSTFdata = theRTVFTobj.visualRFmapPropertiesFromCronerKaplanAnalysis(theLconeVisualRFmap);
    
    % Compute the M-cone visual RF map
    theMconeVisualRFmap = theRTVFTobj.visualRFfromRetinalConePooling(...
                theRTVFTobj.MconeRFcomputeStruct.modelConstants, ...
                theRTVFTobj.MconeRFcomputeStruct.retinalConePoolingParams.finalValues);
  
    % Compute the LMcone visual STF data for the L-cone visual RF map
    theMconeSTFdata = theRTVFTobj.visualRFmapPropertiesFromCronerKaplanAnalysis(theMconeVisualRFmap);
    



    fittedVisualSTF = theLconeSTFdata.fittedDoGModelToVisualSTF.compositeSTF(:);
    theLconeCenterRGCModelVisualSTF = fittedVisualSTF / max(fittedVisualSTF(:));

    fittedVisualSTF = theMconeSTFdata.fittedDoGModelToVisualSTF.compositeSTF(:);
    theMconeCenterRGCModelVisualSTF = fittedVisualSTF / max(fittedVisualSTF(:));


    spatialFrequencySupport = theLconeSTFdata.spatialFrequencySupport(:);
    LconeCenterVisualSTF = theLconeCenterRGCModelVisualSTF * theCenterConeTypeWeights(cMosaic.LCONE_ID)/sum(theCenterConeTypeWeights);
    MconeCenterVisualSTF = theMconeCenterRGCModelVisualSTF * theCenterConeTypeWeights(cMosaic.MCONE_ID)/sum(theCenterConeTypeWeights);
    multiSpectralVisualSTF = LconeCenterVisualSTF + MconeCenterVisualSTF;
end