function generateVisualizationCache(obj, xSupport, ySupport, centerSubregionContourSamples, contourGenerationMethod)

    if (isfield(obj.visualizationCache, 'rfCenterPatchData')) && ...
       (~isempty(obj.visualizationCache.rfCenterPatchData)) && ...
       (isfield(obj.visualizationCache, 'centerSubregionContourSamples')) && ...
       (obj.visualizationCache.centerSubregionContourSamples == centerSubregionContourSamples)
        % Already in visualizationCache, so return
        return;
    end

    obj.visualizationCache = [];
    fprintf('\nComputing RF center outline contours. Please wait ...');
    tic


    spatialSupportSamples = 48;
    if (~isempty(obj.rgcRFcenterConePoolingMatrix))
        minCenterConePoolingWeights = max(obj.rgcRFcenterConePoolingMatrix,[], 1)* 0.001;
    else
        minCenterConePoolingWeights = max(obj.rgcRFcenterConeConnectivityMatrix,[], 1)* 0.001;
    end

    if (~isempty(obj.rgcRFcenterConePoolingMatrix))
        [verticesNumForRGC, verticesList, facesList, colorVertexCData, theContourData, rfCenterConeConnectionLineSegments, singleConeRGCdotPositions] = ...
            graphicDataForSubregion(obj, obj.rgcRFcenterConePoolingMatrix, minCenterConePoolingWeights, ...
            xSupport, ySupport, spatialSupportSamples,centerSubregionContourSamples, contourGenerationMethod);
    else
        [verticesNumForRGC, verticesList, facesList, colorVertexCData, theContourData, rfCenterConeConnectionLineSegments, singleConeRGCdotPositions] = ...
            graphicDataForSubregion(obj, obj.rgcRFcenterConeConnectivityMatrix, minCenterConePoolingWeights, ...
            xSupport, ySupport, spatialSupportSamples,centerSubregionContourSamples, contourGenerationMethod);
    end


    obj.visualizationCache.centerSubregionContourSamples = centerSubregionContourSamples;
    obj.visualizationCache.rfCenterContourData = theContourData;
    obj.visualizationCache.rfCenterPatchData.vertices = verticesList;
    obj.visualizationCache.rfCenterPatchData.verticesNumForRGC = verticesNumForRGC;
    obj.visualizationCache.rfCenterPatchData.faces = facesList;
    obj.visualizationCache.rfCenterPatchData.faceVertexCData = colorVertexCData;
    obj.visualizationCache.rfCenterConeConnectionLineSegments = rfCenterConeConnectionLineSegments;
    obj.visualizationCache.rfCenterSingleConeInputDotPositions = singleConeRGCdotPositions;

    % Find all input cone indices that are connected to the RF centers
    if (~isempty(obj.rgcRFcenterConePoolingMatrix))
        centerConnectedConeIndices = find(sum(obj.rgcRFcenterConePoolingMatrix,2) > 0);
    else
        centerConnectedConeIndices = find(sum(obj.rgcRFcenterConeConnectivityMatrix,2) > 0);
    end

    idx = find(obj.inputConeMosaic.coneTypes(centerConnectedConeIndices) == cMosaic.LCONE_ID);
    obj.visualizationCache.lConeIndicesConnectedToRGCcenters = centerConnectedConeIndices(idx);
    idx = find(obj.inputConeMosaic.coneTypes(centerConnectedConeIndices) == cMosaic.MCONE_ID);
    obj.visualizationCache.mConeIndicesConnectedToRGCcenters= centerConnectedConeIndices(idx);


    if (~isempty(obj.rgcRFsurroundConePoolingMatrix))
        % Find all input cone indices that are connected to the RF surrounds
        surroundConnectedConeIndices = find(sum(obj.rgcRFsurroundConePoolingMatrix,2) > 0);
        xx = squeeze(obj.inputConeMosaic.coneRFpositionsDegs(surroundConnectedConeIndices,1));
        yy = squeeze(obj.inputConeMosaic.coneRFpositionsDegs(surroundConnectedConeIndices,2));
        obj.visualizationCache.surroundConesXrange = [min(xx) max(xx)];
        obj.visualizationCache.surroundConesYrange = [min(yy) max(yy)];
        idx = find(obj.inputConeMosaic.coneTypes(surroundConnectedConeIndices) == cMosaic.LCONE_ID);
        obj.visualizationCache.lConeIndicesConnectedToRGCsurrounds = surroundConnectedConeIndices(idx);
        idx = find(obj.inputConeMosaic.coneTypes(surroundConnectedConeIndices) == cMosaic.MCONE_ID);
        obj.visualizationCache.mConeIndicesConnectedToRGCsurrounds = surroundConnectedConeIndices(idx);
    end

    fprintf(' Done in %2.1f seconds\n', toc);
end



function [verticesNumForRGC, verticesList, facesList, colorVertexCData, theContourData, subregionConeConnectionLineSegments, singleConeRGCdotPositions] = ...
        graphicDataForSubregion(obj, conePoolingMatrix, minPoolingWeights, xSupport, ySupport, spatialSupportSamples, ...
        centerSubregionContourSamples, contourGenerationMethod)
        
    coneApertureSizeSpecifierForRGCRFplotting = 'spacing based';
    %coneApertureSizeSpecifierForRGCRFplotting = 'characteristic radius based';

    switch (coneApertureSizeSpecifierForRGCRFplotting)
        case 'spacing based'
            coneRFradiiDegs = 0.6*0.5*obj.inputConeMosaic.coneRFspacingsDegs;
        case 'characteristic radius based'
            coneRFradiiDegs = ...
                obj.inputConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor * ...
                obj.inputConeMosaic.coneApertureDiametersDegs;
        otherwise
            error('Unknown apertureSizeSpecifierForRGCRFplotting: ''%s''.', coneApertureSizeSpecifierForRGCRFplotting)
    end

    verticesNumForRGC = zeros(1, obj.rgcsNum);
    theContourData = cell(1, obj.rgcsNum);
 
    lineSegmentIndex = 0;
    subregionConeConnectionLineSegments = [];

    singleConeRGCindex = 0;
    singleConeRGCdotPositions = [];

    for iRGC = 1:obj.rgcsNum
        % Retrieve the subregion cone indices & weights
        connectivityVector = full(squeeze(conePoolingMatrix(:, iRGC)));
        subregionConeIndices = find(connectivityVector > minPoolingWeights(iRGC));

        if (numel(subregionConeIndices) == 1)
            singleConeRGCindex = singleConeRGCindex + 1;
            singleConeRGCdotPositions(singleConeRGCindex,:) = obj.inputConeMosaic.coneRFpositionsDegs(subregionConeIndices,:);
        end

        theConePoolingWeights = connectivityVector(subregionConeIndices);
        theConePositions = obj.inputConeMosaic.coneRFpositionsDegs(subregionConeIndices,:);
        theConeRFRadii = coneRFradiiDegs(subregionConeIndices);
        theConeTypes = obj.inputConeMosaic.coneTypes(subregionConeIndices);

        inputConesNum = size(theConePositions,1);
        if (inputConesNum > 1)
            for iCone = 1:inputConesNum
                lineSegmentIndex = lineSegmentIndex + 1;
                subregionConeConnectionLineSegments.Xpos(:, lineSegmentIndex) = ...
                    [obj.rgcRFpositionsDegs(iRGC,1) theConePositions(iCone,1)]';
                subregionConeConnectionLineSegments.Ypos(:, lineSegmentIndex) = ...
                    [obj.rgcRFpositionsDegs(iRGC,2) theConePositions(iCone,2)]';
                subregionConeConnectionLineSegments.coneTypes(lineSegmentIndex) = theConeTypes(iCone);
            end
        end

        switch (contourGenerationMethod)
            case 'ellipseFitBasedOnLocalSpacing'
                rgcRFradius = 0.5*obj.rgcRFspacingsDegs(iRGC);
                theContourData{iRGC} = mRGCMosaic.subregionOutlineContourFromSpacing(...
                    obj.rgcRFpositionsDegs(iRGC,:), rgcRFradius,...
                    xSupport, ySupport, spatialSupportSamples);

            case 'contourOfPooledConeApertureImage'
                theContourData{iRGC} = mRGCMosaic.subregionOutlineContourFromPooledCones(...
                     theConePositions, theConeRFRadii, theConePoolingWeights, ...
                     xSupport, ySupport, spatialSupportSamples);

            case 'ellipseFitToPooledConeApertureImage'
                theContourData{iRGC} = mRGCMosaic.subregionEllipseFromPooledCones(...
                     theConePositions, theConeRFRadii, theConePoolingWeights, ...
                     xSupport, ySupport, spatialSupportSamples, centerSubregionContourSamples);

            otherwise
                error('Unknown contourGenerationMethod: ''%s''.', contourGenerationMethod);
        end % switch

        if (isempty(theContourData{iRGC}))
            verticesNumForRGC(iRGC) = 0;
        else
            s = theContourData{iRGC};
            verticesNumForRGC(iRGC) = size(s.vertices,1);
        end
    end

    maxNoVertices = max(verticesNumForRGC);
    totalVerticesNum = sum(verticesNumForRGC);

    facesList = nan(obj.rgcsNum, maxNoVertices);
    verticesList = zeros(totalVerticesNum,2);
    colorVertexCData = zeros(totalVerticesNum,1);

    currentFacesNum = 0;

    for iRGC = 1:obj.rgcsNum
        if (isempty(theContourData{iRGC}))
            continue;
        end

        s = theContourData{iRGC};
        newVerticesNum = size(s.vertices,1);
        idx = currentFacesNum+(1:newVerticesNum);
        verticesList(idx,:) =  s.vertices;
        colorVertexCData(idx,:) = 0.5;
        facesList(iRGC,1:newVerticesNum) = currentFacesNum + s.faces;
        currentFacesNum = currentFacesNum + newVerticesNum;
    end
end