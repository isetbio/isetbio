function generateVisualizationCache(obj, xSupport, ySupport,  ...
    centerSubregionContourSamples, contourGenerationMethod, ...
    visualizedRGCindices, minConeWeightIncluded,  varargin)

    p = inputParser;
    p.addParameter('forceRegenerateVisualizationCache', false, @islogical);
    p.addParameter('spatialSupportSamples', [], @(x)(isempty(x)||isscalar(x)));
    p.addParameter('maxNumberOfConesOutsideContour', 1, @isscalar);
    p.parse(varargin{:});
    maxNumberOfConesOutsideContour = p.Results.maxNumberOfConesOutsideContour;
    forceRegenerateVisualizationCache = p.Results.forceRegenerateVisualizationCache;
    spatialSupportSamples = p.Results.spatialSupportSamples;
    
    if (isempty(spatialSupportSamples))
        spatialSupportSamples = 48;
    end

    if (isempty(minConeWeightIncluded))
        minConeWeightIncluded = mRGCMosaic.sensitivityAtPointOfOverlap;
    end

    if (isfield(obj.visualizationCache, 'visualizedRGCindices'))
        previousIndicesNum = numel(obj.visualizationCache.visualizedRGCindices);
        visualizedRGCindicesAreIdentical = ...
            ((isempty(previousIndicesNum))&&(isempty(obj.visualizationCache.visualizedRGCindices))) || ...
            ((numel(obj.visualizationCache.visualizedRGCindices)==numel(visualizedRGCindices)))&&(isempty(setdiff(obj.visualizationCache.visualizedRGCindices, visualizedRGCindices)));
    end
    
    conditionsAreTrue = [];
    invalidationReason = {};
    conditionsAreTrue(numel(conditionsAreTrue)+1) = ...
        (isfield(obj.visualizationCache, 'rfCenterPatchData')) && ...
        (~isempty(obj.visualizationCache.rfCenterPatchData));
    invalidationReason{numel(invalidationReason)+1} = sprintf('rfPatchData does not exist or it is empty.');

    conditionsAreTrue(numel(conditionsAreTrue)+1) = ...
        (isfield(obj.visualizationCache, 'centerSubregionContourSamples')) && ...
        (obj.visualizationCache.centerSubregionContourSamples == centerSubregionContourSamples);
    invalidationReason{numel(invalidationReason)+1} = sprintf('centerSubregionsContourSamples does not exist, or a different value is requested.');

    conditionsAreTrue(numel(conditionsAreTrue)+1) = ...
        (isfield(obj.visualizationCache, 'minConeWeightIncluded')) && ...
        (obj.visualizationCache.minConeWeightIncluded == minConeWeightIncluded);
    invalidationReason{numel(invalidationReason)+1} = sprintf('minConeWeightIncluded does not exist, or a different value is requested.');

    conditionsAreTrue(numel(conditionsAreTrue)+1) = ...
        (isfield(obj.visualizationCache, 'maxNumberOfConesOutsideContour')) && ...
        (obj.visualizationCache.maxNumberOfConesOutsideContour == maxNumberOfConesOutsideContour);
    invalidationReason{numel(invalidationReason)+1} = sprintf('maxNumberOfConesOutsideContour does not exist, or a different value is requested.');

    conditionsAreTrue(numel(conditionsAreTrue)+1) = ...
        (isfield(obj.visualizationCache, 'spatialSupportSamples')) && ...
        (obj.visualizationCache.spatialSupportSamples == spatialSupportSamples);
    invalidationReason{numel(invalidationReason)+1} = sprintf('spatialSupportSamples does not exist, or a different value is requested.');


    conditionsAreTrue(numel(conditionsAreTrue)+1) = ...
        (isfield(obj.visualizationCache, 'visualizedRGCindices')) && ...
        (visualizedRGCindicesAreIdentical);
    invalidationReason{numel(invalidationReason)+1} = sprintf('visualizedRGCindices does not exist, or a different value is requested.');

    

    if (all(conditionsAreTrue)) && (~forceRegenerateVisualizationCache)
        return;
    end



    reportReasonCacheIsInvalid = ~true;
    if (reportReasonCacheIsInvalid)
        idx = find(conditionsAreTrue == false);
        for i = 1:numel(idx)
            fprintf('Cache is not valid because %s\n', invalidationReason{idx(i)});
        end
    end

    obj.visualizationCache = [];
    fprintf('\nComputing RF center outline contours. Please wait ...');
    tic

    [verticesNumForRGC, verticesList, facesList, colorVertexCData, ...
     theContourData, rfCenterConeConnectionLineSegments, singleConeRGCdotPositions] = graphicDataForSubregion(obj, ...
        obj.rgcRFcenterConeConnectivityMatrix, minConeWeightIncluded, maxNumberOfConesOutsideContour, ...
        xSupport, ySupport, spatialSupportSamples, centerSubregionContourSamples, ...
        contourGenerationMethod, visualizedRGCindices);
    

    obj.visualizationCache.visualizedRGCindices = visualizedRGCindices;
    obj.visualizationCache.centerSubregionContourSamples = centerSubregionContourSamples;
    obj.visualizationCache.minConeWeightIncluded = minConeWeightIncluded;
    obj.visualizationCache.maxNumberOfConesOutsideContour = maxNumberOfConesOutsideContour;
    obj.visualizationCache.spatialSupportSamples = spatialSupportSamples;
    obj.visualizationCache.rfCenterContourData = theContourData;
    obj.visualizationCache.rfCenterPatchData.vertices = verticesList;
    obj.visualizationCache.rfCenterPatchData.verticesNumForRGC = verticesNumForRGC;
    obj.visualizationCache.rfCenterPatchData.faces = facesList;
    obj.visualizationCache.rfCenterPatchData.faceVertexCData = colorVertexCData;
    obj.visualizationCache.rfCenterConeConnectionLineSegments = rfCenterConeConnectionLineSegments;
    obj.visualizationCache.rfCenterSingleConeInputDotPositions = singleConeRGCdotPositions;


    if (1==2)
        % NOT SURE IF WE NEED THE BELOW, so (1==2)

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

    end

    fprintf(' Done in %2.1f seconds\n', toc);
end


function [verticesNumForRGC, verticesList, facesList, colorVertexCData, ...
         theContourData, subregionConeConnectionLineSegments, singleConeRGCdotPositions] = graphicDataForSubregion(obj, ...
            conePoolingMatrix, minConeWeightIncluded, maxNumberOfConesOutsideContour, xSupport, ySupport, spatialSupportSamples, ...
            centerSubregionContourSamples, contourGenerationMethod, visualizedRGCindices)
        
    verticesNumForRGC = zeros(1, numel(visualizedRGCindices));
    theContourData = cell(1, numel(visualizedRGCindices));

    lineSegmentIndex = 0;
    subregionConeConnectionLineSegments = [];

    singleConeRGCindex = 0;
    singleConeRGCdotPositions = [];

    skipConnectivityVectorComputation = strcmp(contourGenerationMethod, 'ellipseFitBasedOnLocalSpacingOfSourceLattice');

    parfor iRGC = 1:numel(visualizedRGCindices)
        theRGCindex = visualizedRGCindices(iRGC);

        if (skipConnectivityVectorComputation)
            rgcRFradius = 0.5*obj.rgcRFspacingsDegsOfSourceLattice(theRGCindex);
            theContourData{iRGC} = subregionOutlineContourFromSpacing(...
                        obj.rgcRFpositionsDegsOfSourceLattice(theRGCindex,:), rgcRFradius,...
                        xSupport, ySupport, spatialSupportSamples);
            continue;
        end

        % Retrieve the subregion cone indices & weights
        connectivityVector = full(squeeze(conePoolingMatrix(:, theRGCindex)));
        subregionConeIndices = find(connectivityVector >  minConeWeightIncluded);

        if (numel(subregionConeIndices) > 0)
            % For visualization purposes, all included cones have unit amplitude
            theConePoolingWeights = 1 + 0*connectivityVector(subregionConeIndices);
            theConePositions = obj.inputConeMosaic.coneRFpositionsDegs(subregionConeIndices,:);
            theConeSpacings = obj.inputConeMosaic.coneRFspacingsDegs(subregionConeIndices);

            switch (contourGenerationMethod)
                case 'ellipseFitBasedOnLocalSpacing'
                    rgcRFradius = 0.5*obj.rgcRFspacingsDegs(theRGCindex);
                    theContourData{iRGC} = subregionOutlineContourFromSpacing(...
                        obj.rgcRFpositionsDegs(theRGCindex,:), rgcRFradius,...
                        xSupport, ySupport, spatialSupportSamples);
    
                case 'contourOfPooledConeApertureImage'
                    theContourData{iRGC} = subregionOutlineContourFromPooledCones(...
                         theConePositions, theConeSpacings, theConePoolingWeights, ...
                         xSupport, ySupport, spatialSupportSamples);
    
                case 'ellipseFitToPooledConeApertureImage'
                    theContourData{iRGC} = mRGCMosaic.subregionEllipseFromPooledCones(...
                         theConePositions, theConeSpacings, theConePoolingWeights, ...
                         xSupport, ySupport, spatialSupportSamples, centerSubregionContourSamples);
                    if (isempty(theContourData{iRGC}))
                        fprintf(2, 'poolingWeights for RGC #%d is []. Returning an empty contourData struct\n', theRGCindex);
                    end
    
                case 'ellipseFitToPooledConePositions'
                    inputConeNumerosityForThisRGC = size(theConePositions,1);
                    minInputConeNumerosityToEmployEllipseFitToPooledConePositions = 5;
                    if (inputConeNumerosityForThisRGC >= minInputConeNumerosityToEmployEllipseFitToPooledConePositions)
                        ellipseContourAngles = 0:10:350;
                        % pLevel = 0.8; % The higher, the wider the Gaussian
                        % Or [], and specify # of cones to leave outside of the ellipse
                        pLevel = [];  % When empty, return a contour that does not include the maxNumberOfConesOutsideContour furthest cones
                        theContourData{iRGC} = mRGCMosaic.subregionEllipseFromPooledConePositions(...
                             theConePositions, centerSubregionContourSamples, ellipseContourAngles, pLevel, maxNumberOfConesOutsideContour);
                    else
                        
                        theContourData{iRGC} = mRGCMosaic.subregionEllipseFromPooledCones(...
                         theConePositions, theConeSpacings, theConePoolingWeights, ...
                         xSupport, ySupport, spatialSupportSamples, centerSubregionContourSamples);
                    end

                    if (isempty(theContourData{iRGC}))
                        fprintf(2, 'poolingWeights for RGC #%d is []. Returning an empty contourData struct\n', theRGCindex);
                    end

                otherwise
                    error('Unknown contourGenerationMethod: ''%s''.', contourGenerationMethod);
            end % switch
        else
            fprintf('RGC #%d contains no cones with a weight > min weight specified (%f) and will therefore not be visualized. \n', theRGCindex, minConeWeightIncluded);
            theContourData{iRGC} = [];
        end

        if (isempty(theContourData{iRGC}))
            verticesNumForRGC(iRGC) = 0;
        else
            s = theContourData{iRGC};
            verticesNumForRGC(iRGC) = size(s.vertices,1);
        end
    end % parfor iRGC

    for iRGC = 1:numel(visualizedRGCindices)
        theRGCindex = visualizedRGCindices(iRGC);
        
        if (skipConnectivityVectorComputation)
            continue;
        end

        % Retrieve the subregion cone indices & weights
        connectivityVector = full(squeeze(conePoolingMatrix(:, theRGCindex)));
        subregionConeIndices = find(connectivityVector >  minConeWeightIncluded);
        
        if (subregionConeIndices > 0)
            % Compute the subregion centroid (weighted by the connectivity vector)
            subregionCentroidDegs = sum(connectivityVector(subregionConeIndices) .* obj.inputConeMosaic.coneRFpositionsDegs(subregionConeIndices,:),1) / ...
                                    sum(connectivityVector(subregionConeIndices));

            if (numel(subregionConeIndices) == 1)
                singleConeRGCindex = singleConeRGCindex + 1;
                singleConeRGCdotPositions(singleConeRGCindex,:) = obj.inputConeMosaic.coneRFpositionsDegs(subregionConeIndices,:);
            end

            theConePositions = obj.inputConeMosaic.coneRFpositionsDegs(subregionConeIndices,:);
            theConeTypes = obj.inputConeMosaic.coneTypes(subregionConeIndices);
            theConePoolingWeights = connectivityVector(subregionConeIndices);

            inputConesNum = size(theConePositions,1);
            if (inputConesNum > 1)
                for iCone = 1:inputConesNum
                    lineSegmentIndex = lineSegmentIndex + 1;
                    subregionConeConnectionLineSegments.Xpos(:, lineSegmentIndex) = ...
                        [subregionCentroidDegs(1) theConePositions(iCone,1)]';
                    subregionConeConnectionLineSegments.Ypos(:, lineSegmentIndex) = ...
                        [subregionCentroidDegs(2) theConePositions(iCone,2)]';
                    subregionConeConnectionLineSegments.lineSegmentWidths(lineSegmentIndex) = theConePoolingWeights(iCone)/max(theConePoolingWeights(:)) * 4;
                    subregionConeConnectionLineSegments.coneTypes(lineSegmentIndex) = theConeTypes(iCone);
                end
            end 
        end
    end % for iRGC

    maxNoVertices = max(verticesNumForRGC);
    totalVerticesNum = sum(verticesNumForRGC);

    facesList = nan(numel(visualizedRGCindices), maxNoVertices);
    verticesList = zeros(totalVerticesNum,2);
    colorVertexCData = zeros(totalVerticesNum,1);

    currentFacesNum = 0;

    for iRGC = 1:numel(visualizedRGCindices)
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
        xSupport = linspace(min(xx)-10*xSep,max(xx)+10*xSep,spatialSupportSamples);
    end

    if (isempty(ySupport))
        yy = rgcPos(:,2);
        ySupport = linspace(min(yy)-10*xSep,max(yy)+10*xSep,spatialSupportSamples);
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
    s = contourDataFrom2DRFmap(spatialSupportXY, RF, zLevels);
    contourData = s{1};
end

function [contourData, RF] = subregionOutlineContourFromPooledCones(...
    conePos, coneSpacings, poolingWeights, ...
    xSupport, ySupport, spatialSupportSamples, centerSubregionContourSamples)

    if (isempty(poolingWeights))
        contourData = [];
        return;
    end

    % Compute spatial support
    xSep = max(coneSpacings)*sqrt(numel(poolingWeights));
    if (isempty(xSupport))
        xx = conePos(:,1);
        xSupport = linspace(min(xx)-10*xSep,max(xx)+10*xSep,spatialSupportSamples);
    end

    if (isempty(ySupport))
        yy = conePos(:,2);
        ySupport = linspace(min(yy)-10*xSep,max(yy)+10*xSep,spatialSupportSamples);
    end

    [X,Y] = meshgrid(xSupport, ySupport);
    spatialSupportXY(:,1) = xSupport(:);
    spatialSupportXY(:,2) = ySupport(:);

    RF = zeros(size(X));
    for iCone = 1:numel(poolingWeights)
        % Characteristic radius of the input RF
        rC = exp(-1) * coneSpacings(iCone);
        % Compute aperture2D x weight
        XX = X-conePos(iCone,1);
        YY = Y-conePos(iCone,2);
        theAperture2D = poolingWeights(iCone) * exp(-(XX/rC).^2) .* exp(-(YY/rC).^2);
        % Accumulate 2D apertures
        RF = RF + theAperture2D;
    end

    zLevels(1) = 0.05*min(poolingWeights);
    zLevels(2) = 0.06*min(poolingWeights);
    RF(RF<zLevels(1)) = 0.0;
    RF(RF>=zLevels(1)) = 0.1;
    RF(RF>=zLevels(2)) = 0.11;
    zLevels = [0.1 0.11];
    s = contourDataFrom2DRFmap(spatialSupportXY, RF, zLevels);
    contourData = s{1};
end



function cData = contourDataFrom2DRFmap(spatialSupportXY, zData, zLevels)
    xSupport = squeeze(spatialSupportXY(:,1));
    ySupport = squeeze(spatialSupportXY(:,2));
    C = contourc(xSupport, ySupport, zData, zLevels);
    dataPoints = size(C,2);
    startPoint = 1;

    iContour = 0;
    levelsNum = 0;
    while (startPoint < dataPoints)
        levelsNum = levelsNum + 1;
        if (levelsNum > 1)
            theLevel = C(1,startPoint);
        end
        theLevelVerticesNum = C(2,startPoint);
        x = C(1,startPoint+(1:theLevelVerticesNum));
        y = C(2,startPoint+(1:theLevelVerticesNum));
        v = [x(:) y(:)];

        f = 1:numel(x);
        iContour = iContour+1;
        cData{iContour} = struct('faces', f, 'vertices', v);
        startPoint = startPoint + theLevelVerticesNum+1;
    end
end


