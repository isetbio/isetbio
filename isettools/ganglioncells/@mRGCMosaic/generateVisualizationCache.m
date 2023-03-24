function generateVisualizationCache(obj, xSupport, ySupport)
    if (isfield(obj.visualizationCache, 'rfCenterPatchData')) && ...
       (~isempty(obj.visualizationCache.rfCenterPatchData))
        % Already in visualizationCache, so return
        return;
    end

    fprintf('\nComputing RF center outline contours. Please wait ...');
    tic

    % Compute graphic data for center contours
    spatialSupportSamples = 24;
    
    if (~isempty(obj.rgcRFcenterConePoolingMatrix))
        minCenterConePoolingWeights = max(obj.rgcRFcenterConePoolingMatrix,[], 1)* 0.001;
    else
        minCenterConePoolingWeights = max(obj.rgcRFcenterConeConnectivityMatrix,[], 1)* 0.001;
    end

    if (~isempty(obj.rgcRFcenterConePoolingMatrix))
        [verticesNumForRGC, verticesList, facesList, colorVertexCData, theContourData, rfCenterConeConnectionLineSegments] = ...
            graphicDataForSubregion(obj, obj.rgcRFcenterConePoolingMatrix, minCenterConePoolingWeights, ...
            xSupport, ySupport, spatialSupportSamples);
    else
        [verticesNumForRGC, verticesList, facesList, colorVertexCData, theContourData, rfCenterConeConnectionLineSegments] = ...
            graphicDataForSubregion(obj, obj.rgcRFcenterConeConnectivityMatrix, minCenterConePoolingWeights, ...
            xSupport, ySupport, spatialSupportSamples);
    end



    obj.visualizationCache.rfCenterContourData = theContourData;
    obj.visualizationCache.rfCenterPatchData.vertices = verticesList;
    obj.visualizationCache.rfCenterPatchData.verticesNumForRGC = verticesNumForRGC;
    obj.visualizationCache.rfCenterPatchData.faces = facesList;
    obj.visualizationCache.rfCenterPatchData.faceVertexCData = colorVertexCData;
    obj.visualizationCache.rfCenterConeConnectionLineSegments = rfCenterConeConnectionLineSegments;

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

function [verticesNumForRGC, verticesList, facesList, colorVertexCData, theContourData, subregionConeConnectionLineSegments] = ...
        graphicDataForSubregion(obj, conePoolingMatrix, minPoolingWeights, xSupport, ySupport, spatialSupportSamples)
        
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
    for iRGC = 1:obj.rgcsNum
        % Retrieve the subregion cone indices & weights
        connectivityVector = full(squeeze(conePoolingMatrix(:, iRGC)));
        subregionConeIndices = find(connectivityVector > minPoolingWeights(iRGC));

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


        theContourData{iRGC} = subregionOutlineContourFromPooledCones(...
                 theConePositions, theConeRFRadii, theConePoolingWeights, ...
                 xSupport, ySupport, spatialSupportSamples);

        s = theContourData{iRGC}{1};
        verticesNumForRGC(iRGC) = size(s.vertices,1);
    end

    maxNoVertices = max(verticesNumForRGC);
    totalVerticesNum = sum(verticesNumForRGC);

    facesList = nan(obj.rgcsNum, maxNoVertices);
    verticesList = zeros(totalVerticesNum,2);
    colorVertexCData = zeros(totalVerticesNum,1);

    currentFacesNum = 0;

    for iRGC = 1:obj.rgcsNum
        s = theContourData{iRGC}{1};
        newVerticesNum = size(s.vertices,1);
        idx = currentFacesNum+(1:newVerticesNum);
        verticesList(idx,:) =  s.vertices;
        colorVertexCData(idx,:) = 0.5;
        facesList(iRGC,1:newVerticesNum) = currentFacesNum + s.faces;
        currentFacesNum = currentFacesNum + newVerticesNum;
    end

end


function contourData = subregionOutlineContourFromPooledCones(...
    conePos, coneRc, poolingWeights, ...
    xSupport, ySupport, spatialSupportSamples)

    % Compute spatial support
    xSep = max(coneRc)*2*sqrt(numel(poolingWeights));
    if (isempty(xSupport))
        xx = conePos(:,1);
        xSupport = linspace(min(xx)-xSep,max(xx)+xSep,spatialSupportSamples);
    end

    if (isempty(ySupport))
        yy = conePos(:,2);
        ySupport = linspace(min(yy)-xSep,max(yy)+xSep,spatialSupportSamples);
    end

    [X,Y] = meshgrid(xSupport, ySupport);
    spatialSupportXY(:,1) = xSupport(:);
    spatialSupportXY(:,2) = ySupport(:);

    RF = zeros(size(X));
    for iCone = 1:numel(poolingWeights)
        % Characteristic radius of the input RF
        rC = coneRc(iCone);
        % Compute aperture2D x weight
        XX = X-conePos(iCone,1);
        YY = Y-conePos(iCone,2);
        theAperture2D = poolingWeights(iCone) * exp(-(XX/rC).^2) .* exp(-(YY/rC).^2);
        % Accumulate 2D apertures
        RF = RF + theAperture2D;
    end

    zLevels(1) = 0.02*min(poolingWeights);
    zLevels(2) = max(poolingWeights);

    contourData = mRGCMosaic.contourDataFromDensityMap(spatialSupportXY, RF, zLevels);
end



