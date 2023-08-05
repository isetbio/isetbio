function generateVisualizationCache(obj, xSupport, ySupport, centerSubregionContourSamples, contourGenerationMethod)

    if (isfield(obj.visualizationCache, 'rfCenterPatchData')) && ...
       (~isempty(obj.visualizationCache.rfCenterPatchData)) && ...
       (isfield(obj.visualizationCache, 'obj.visualizationCache')) && ...
       (obj.visualizationCache.centerSubregionContourSamples == centerSubregionContourSamples)
        % Already in visualizationCache, so return
        return;
    end

    fprintf('\nComputing RF center outline contours. Please wait ...');
    tic


    spatialSupportSamples = 48;
    if (~isempty(obj.rgcRFcenterConePoolingMatrix))
        minCenterConePoolingWeights = max(obj.rgcRFcenterConePoolingMatrix,[], 1)* 0.001;
    else
        minCenterConePoolingWeights = max(obj.rgcRFcenterConeConnectivityMatrix,[], 1)* 0.001;
    end

    if (~isempty(obj.rgcRFcenterConePoolingMatrix))
        [verticesNumForRGC, verticesList, facesList, colorVertexCData, theContourData, rfCenterConeConnectionLineSegments] = ...
            graphicDataForSubregion(obj, obj.rgcRFcenterConePoolingMatrix, minCenterConePoolingWeights, ...
            xSupport, ySupport, spatialSupportSamples,centerSubregionContourSamples, contourGenerationMethod);
    else
        [verticesNumForRGC, verticesList, facesList, colorVertexCData, theContourData, rfCenterConeConnectionLineSegments] = ...
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

        switch (contourGenerationMethod)
            case 'ellipseFitBasedOnLocalSpacing'
                rgcRFradius = 0.5*obj.rgcRFspacingsDegs(iRGC);
                theContourData{iRGC} = subregionOutlineContourFromSpacing(...
                    obj.rgcRFpositionsDegs(iRGC,:), rgcRFradius,...
                    xSupport, ySupport, spatialSupportSamples);
            case 'contourOfPooledConeApertureImage'
                theContourData{iRGC} = subregionOutlineContourFromPooledCones(...
                     theConePositions, theConeRFRadii, theConePoolingWeights, ...
                     xSupport, ySupport, spatialSupportSamples);

            case 'ellipseFitToPooledConeApertureImage'
                theContourData{iRGC} = subregionEllipseFromPooledCones(...
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


function contourData = subregionOutlineContourFromSpacing(...
    rgcPos, rgcRFradius,  ...
    xSupport, ySupport, spatialSupportSamples)

    % Compute spatial support
    xSep = rgcRFradius;
    if (isnan(xSep))
        contourData = [];
        return;
    end

    if (isempty(xSupport))
        xx = rgcPos(:,1);
        xSupport = linspace(min(xx)-xSep,max(xx)+xSep,spatialSupportSamples);
    end

    if (isempty(ySupport))
        yy = rgcPos(:,2);
        ySupport = linspace(min(yy)-xSep,max(yy)+xSep,spatialSupportSamples);
    end

    [X,Y] = meshgrid(xSupport, ySupport);
    spatialSupportXY(:,1) = xSupport(:);
    spatialSupportXY(:,2) = ySupport(:);
    
    % 2*sigma = radius
    rSigma = rgcRFradius/3;
    % Compute aperture2D x weight
    XX = X-rgcPos(1);
    YY = Y-rgcPos(2);
    theAperture2D = exp(-0.5*(XX/rSigma).^2) .* exp(-0.5*(YY/rSigma).^2);
    RF = theAperture2D;

    zLevels(1) = exp(-2.5);
    zLevels(2) = exp(-2.4);
    s = mRGCMosaic.contourDataFromDensityMap(spatialSupportXY, RF, zLevels);
    contourData = s{1};
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
    s = mRGCMosaic.contourDataFromDensityMap(spatialSupportXY, RF, zLevels);
    contourData = s{1};
end

function contourData = subregionEllipseFromPooledCones(...
    conePos, coneRc, poolingWeights, ...
    xSupport, ySupport, spatialSupportSamples, centerSubregionContourSamples)

    % Compute spatial support
    xSep = max(coneRc)*2.5*sqrt(numel(poolingWeights));

    if (isempty(poolingWeights))
        fprintf(2, 'poolingWeights is []. Returning an empty contourData struct\n');
        contourData = [];
        return;
    end

    if (isempty(xSupport))
        xx = conePos(:,1);
        xSupport = linspace(min(xx)-xSep,max(xx)+xSep,spatialSupportSamples);
    end

    if (isempty(ySupport))
        yy = conePos(:,2);
        ySupport = linspace(min(yy)-xSep,max(yy)+xSep,spatialSupportSamples);
    end

    [X,Y] = meshgrid(xSupport, ySupport);
    RF = zeros(size(X));

    if (numel(poolingWeights)>1)
        xx = conePos(:,1);
        yy = conePos(:,2);
        dx = max(xx)-min(xx);
        dy = max(yy)-min(yy);
        if (dx > dy)
            [~,idx] = sort(xx);
        else
            [~,idx] = sort(yy);
        end
        
        xo = mean(xx);
        yo = mean(yy);

        poolingWeights = poolingWeights(idx);
        coneRc = coneRc(idx);
        xx = xx(idx);
        yy = yy(idx);
        interpolationPointsNum = 4;
        xxxInterp = zeros(1, interpolationPointsNum * numel(xx));
        yyyInterp = zeros(1, interpolationPointsNum * numel(xx));
        coneRcInterp = zeros(1, interpolationPointsNum * numel(xx));
        poolingWeightsInterp = zeros(1, interpolationPointsNum * numel(xx));

        for ii = 1:numel(xx)
            for iii = 1:interpolationPointsNum
                f = (iii-1)/interpolationPointsNum;
                xxxInterp((ii-1)*interpolationPointsNum+iii) = xo*(1-f) + xx(ii)*f;
                yyyInterp((ii-1)*interpolationPointsNum+iii) = yo*(1-f) + yy(ii)*f;
                coneRcInterp((ii-1)*interpolationPointsNum+iii) = coneRc(ii);
                poolingWeightsInterp((ii-1)*interpolationPointsNum+iii) = poolingWeights(ii);
            end
        end


        for iConeInterp = 1:numel(coneRcInterp)
            % Characteristic radius of the input RF
            rC = coneRcInterp(iConeInterp);
            % Compute aperture2D x weight
            XX = X-xxxInterp(iConeInterp);
            YY = Y-yyyInterp(iConeInterp);
            theAperture2D = poolingWeightsInterp(iConeInterp) * exp(-(XX/rC).^2) .* exp(-(YY/rC).^2);
            % Accumulate 2D apertures
            RF = RF + theAperture2D;
            RF(RF>1) = 1;
        end

        % for iCone = 1:numel(poolingWeights)
        %     % Characteristic radius of the input RF
        %     rC = coneRc(iCone);
        %     % Compute aperture2D x weight
        %     XX = X-conePos(iCone,1);
        %     YY = Y-conePos(iCone,2);
        %     theAperture2D = poolingWeights(iCone) * exp(-(XX/rC).^2) .* exp(-(YY/rC).^2);
        %     % Accumulate 2D apertures
        %     RF = RF + theAperture2D;
        % end
    else
        % Characteristic radius of the input RF
        iCone = 1;
        rC = coneRc(iCone);
        % Compute aperture2D x weight
        XX = X-conePos(iCone,1);
        YY = Y-conePos(iCone,2);
        theAperture2D = poolingWeights(iCone) * exp(-(XX/rC).^2) .* exp(-(YY/rC).^2);
        RF = theAperture2D;
        RF(RF>1) = 1;
    end


    % Binarize
    RF = RF / max(RF(:));
    RF(RF<0.1) = 0.0;
    RF(RF>0) = 1.0;
    BW = imbinarize(RF);


    % Extract the maximum area
    BW = imclearborder(BW);
    BW = bwareafilt(BW,1);

    % Calculate centroid, orientation and major/minor axis length of the ellipse
    s = regionprops(BW,{'Centroid','Orientation','MajorAxisLength','MinorAxisLength'});
    if (isempty(s))
       
        figure()
        subplot(1,2,1);
        imagesc(RF)
        axis 'image'
        subplot(1,2,2)
        imagesc(BW);
        axis 'image'
        pause
    end

    % Calculate the ellipse line
    theta = linspace(0, 2*pi, centerSubregionContourSamples);
    col = (s.MajorAxisLength/2)*cos(theta);
    row = (s.MinorAxisLength/2)*sin(theta);
    M = makehgtform('translate',[s.Centroid, 0],'zrotate',deg2rad(-1*s.Orientation));
    D = M*[col;row;zeros(1,numel(row));ones(1,numel(row))];

    x = D(1,:);
    y = D(2,:);
    x = (x-1)/(numel(xSupport)-1) * (xSupport(end)-xSupport(1)) + xSupport(1); 
    y = (y-1)/(numel(ySupport)-1) * (ySupport(end)-ySupport(1)) + ySupport(1); 

    v = [x(:) y(:)];
    f = 1:numel(x);
    s = struct('faces', f, 'vertices', v);
    contourData = s;
end






