function visualizeVisualRFcentersOfComputeReadyMidgetRGCMosaic(...
    theComputeReadyMRGCmosaic, optimallyMappedSubspaceRFmapsFileName)

    % Load tha optimally mapped visual RF maps.
    % These are computed via
    % MosaicPoolingOptimizer.computeVisualRFsOfComputeReadyMidgetRGCMosaic()
    load(optimallyMappedSubspaceRFmapsFileName, 'optimallyMappedVisualRFmaps');

    % Generate contour data for all mapped RFs
    zLevels = exp(-1)*[0.99 1];
    verticesNumForRGC = zeros(1, numel(optimallyMappedVisualRFmaps));
    theContourData = cell(1, numel(optimallyMappedVisualRFmaps));

    for iRGC = 1:numel(optimallyMappedVisualRFmaps)
        s = optimallyMappedVisualRFmaps{iRGC};
        if (isempty(s))
            % optimally mapped RF not computed yet
            continue;
        end

        spatialSupportXY = [s.spatialSupportDegsX(:) s.spatialSupportDegsY(:)];
        theContourData{iRGC} = mRGCMosaic.contourDataFromDensityMap(spatialSupportXY, double(s.theRFmap), zLevels);

        s = theContourData{iRGC}{1};
        verticesNumForRGC(iRGC) = size(s.vertices,1);
    end

    maxNoVertices = max(verticesNumForRGC);
    totalVerticesNum = sum(verticesNumForRGC);

    facesList = nan(numel(optimallyMappedVisualRFmaps), maxNoVertices);
    verticesList = zeros(totalVerticesNum,2);
    colorVertexCData = zeros(totalVerticesNum,1);

    currentFacesNum = 0;

    for iRGC = 1:numel(optimallyMappedVisualRFmaps)
        if (~isempty(theContourData{iRGC}))
            s = theContourData{iRGC}{1};
            newVerticesNum = size(s.vertices,1);
            idx = currentFacesNum+(1:newVerticesNum);
            verticesList(idx,:) =  s.vertices;
            colorVertexCData(idx,:) = 0.5;
            facesList(iRGC,1:newVerticesNum) = currentFacesNum + s.faces;
            currentFacesNum = currentFacesNum + newVerticesNum;
        end
    end

    rfCenterPatchData.vertices = verticesList;
    rfCenterPatchData.verticesNumForRGC = verticesNumForRGC;
    rfCenterPatchData.faces = facesList;
    rfCenterPatchData.faceVertexCData = colorVertexCData;

    hFig = figure(100); clf;
    set(hFig, 'Position', [10 10 1800 1000], 'Color', [1 1 1]);
    ax = subplot('Position', [0.02 0.02 0.98 0.98]);

    S.Vertices = rfCenterPatchData.vertices;
    S.Faces = rfCenterPatchData.faces;
    S.FaceVertexCData = rfCenterPatchData.faceVertexCData;
    
    S.FaceColor = 'flat';
    S.FaceColor = [0.95 0.95 0.95]*0.6;
    S.EdgeColor = [0 0 0];
        
    S.FaceAlpha = 0.4;
    S.EdgeAlpha = 0.0;
    S.LineWidth = 2;
    patch(S, 'Parent', ax)


%     if (~isempty(labelRGCsWithIndices))
%        hold(ax, 'on')
%        for iRGC = 1:numel(labelRGCsWithIndices)
%                 theRGCindex = labelRGCsWithIndices(iRGC);
%                 S = obj.visualizationCache.rfCenterContourData{theRGCindex}{1};
%                 S.FaceVertexCData = 0.5;
%                 S.FaceColor = 'flat';
%                 S.EdgeColor = [1 1 0];
%                 S.FaceAlpha = 0.0;
%                 S.LineWidth = 2.0;
%                 S.LineStyle = '-';
%                 patch(S, 'Parent', ax);
%        end
%     end

end