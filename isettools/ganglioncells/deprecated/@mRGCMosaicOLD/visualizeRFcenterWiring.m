% Static method to visualize the RF center wiring (fast)
% It requires the connections [N x 2] matrix which contains all the
% pairs [RGCindex coneIndex] of cone indices that feed into each RGC index
        
function hFig = visualizeRFcenterWiring(RGCpos, RGCspacings, conePos, coneSpacings, coneTypes, connections, targetEcc)
    
    hFig = figure(1); clf;
    axesHandle = axes('Position',[0.05 0.05 0.9 0.9]);
    
    cMap = [1 0 0; 0 1 0; 0 0 1; 0.5 0.5 0.5];
    colormap(axesHandle, cMap);
    hold(axesHandle, 'on');
    
    % Find RGC closest to target Ecc
    d = bsxfun(@minus, RGCpos, targetEcc);
    [~, targetRGC] = min(sum(d.^2,2));
    visualizeSingleRGCWiring(axesHandle, RGCpos, RGCspacings, conePos, coneSpacings, coneTypes, connections, targetRGC);

    
    d = sqrt(sum((bsxfun(@minus, conePos, RGCpos(targetRGC,:))).^2,2));
    dxdy = abs(bsxfun(@minus, conePos, RGCpos(targetRGC,:)));
    RGCecc = sqrt(sum(RGCpos(targetRGC,:).^2,2));
    xyRange = max([0.1 0.1+0.05*RGCecc]);
    
    % Plot nearby cones (whether connected or not)
    %nearbyConeIndices = find(d <= xyRange*1.1);
    nearbyConeIndices = find((dxdy(:,1)< xyRange*1.1) & (dxdy(:,2)< xyRange*1.1));
    
    conesAreInputsToSingleRGC = false;
    plotApertures(axesHandle, [], [], conePos, coneSpacings, coneTypes, [], nearbyConeIndices, 1.0, 0.5, conesAreInputsToSingleRGC);
    

    % Plot cones connected to targetRGCs
    plottedRGCs = 'allInNeighborhood';  % Choose from {'allInNeighborhood', 'scattered'};
    
    switch plottedRGCs
        case 'allInNeighborhood'
            d = sqrt(sum((bsxfun(@minus, RGCpos, targetEcc)).^2,2));
            dxdy = abs(bsxfun(@minus, RGCpos, targetEcc));
            targetRGCs = find((dxdy(:,1) < xyRange)  & (dxdy(:,2) < xyRange));
            for k = 1:numel(targetRGCs)
                visualizeSingleRGCWiring(axesHandle, RGCpos, RGCspacings, conePos, coneSpacings, coneTypes, connections, targetRGCs(k));
            end
        case 'scattered'
            for dx = -xyRange/2 : xyRange*0.5 : xyRange/2
                for dy = -xyRange/2 : xyRange*0.5 : xyRange/2
                    d = bsxfun(@minus, RGCpos, targetEcc + [dx dy]);
                    [~, targetRGC] = min(sum(d.^2,2));
                    visualizeSingleRGCWiring(axesHandle, RGCpos, RGCspacings, conePos, coneSpacings, coneTypes, connections, targetRGC);
                end
            end
    end   
end

function visualizeSingleRGCWiring(axesHandle, RGCpos, RGCspacings, conePos, coneSpacings, coneTypes, connections, targetRGC)
   
    % Find cones connected to this RGC RF center
    rgcIDs = connections(:,1);
    coneIDs = connections(:,2);
    
    idx = find(rgcIDs == targetRGC);
    if (isempty(idx))
        return;
    end
    
    connectedConeIndices = coneIDs(idx);
    plotApertures(axesHandle, RGCpos, RGCspacings, conePos, coneSpacings, coneTypes, targetRGC, connectedConeIndices, 1.0, 0.5, true);
    
end

function plotApertures(axesHandle, RGCpos, RGCspacings, conePos, coneSpacings, coneTypes, targetRGC, coneIndices, lineWidth, faceAlpha, conesAreInputsToSingleRGC)
    
    % Retrieve cone apertures and positions
    coneApertureRadii = 0.8*0.5*coneSpacings(coneIndices);
    coneAperturePositions = conePos(coneIndices,:);
    
    if (conesAreInputsToSingleRGC)
        xyMin = min(coneAperturePositions,[],1);
        xyMax = max(coneAperturePositions,[],1);

        m = 4*mean(coneApertureRadii);
        xSupport = linspace(xyMin(1)-m, xyMax(1)+m, 30);
        ySupport = linspace(xyMin(2)-m, xyMax(2)+m, 30);
        [X,Y] = meshgrid(xSupport,ySupport);
        sigma = m/4;
        p = 2;
        for k = 1:numel(coneApertureRadii)
            xo = coneAperturePositions(k,1);
            yo = coneAperturePositions(k,2);
            z = exp(-0.5*((X-xo)/sigma).^p) .* exp(-0.5*((Y-yo)/sigma).^p);
            if (k == 1)
                rf = z;
            else
                rf = rf + z;
            end
        end
        rf = rf / max(rf(:));
   
    else
        rf = [];
    end
    
    % For the L-cones
    lConeIndices = find(coneTypes(coneIndices) == cMosaic.LCONE_ID);
    lConeApertureRadii = coneApertureRadii(lConeIndices);
    lConeAperturePositions = coneAperturePositions(lConeIndices,:);
    
    % For the M-cones
    mConeIndices = find(coneTypes(coneIndices) == cMosaic.MCONE_ID);
    mConeApertureRadii = coneApertureRadii(mConeIndices);
    mConeAperturePositions = coneAperturePositions(mConeIndices,:);
    
    % For the S-cones
    sConeIndices = find(coneTypes(coneIndices) == cMosaic.SCONE_ID);
    sConeApertureRadii = coneApertureRadii(sConeIndices);
    sConeAperturePositions = coneAperturePositions(sConeIndices,:);
    
    
    % Plot
    deltaAngle = 15;
    iTheta = (0:deltaAngle:360) / 180 * pi;
    coneApertureShape.x = cos(iTheta);
    coneApertureShape.y = sin(iTheta);

    % Plot params
    edgeColor = 'none'; %[0 0 0];  
    lineStyle = '-';
     
    
    % Plot L-cones
    if (~isempty(lConeIndices))
        renderPatchArray(axesHandle, coneApertureShape, lConeApertureRadii, ...
            lConeAperturePositions, 1/4*0.9, edgeColor, lineWidth, lineStyle, faceAlpha);
    end
    
    % Plot M-cones
    if (~isempty(mConeIndices))
        renderPatchArray(axesHandle, coneApertureShape, mConeApertureRadii, ...
            mConeAperturePositions, 2/4*0.9, edgeColor, lineWidth, lineStyle, faceAlpha);
    end
    
    % Plot S-cones
    if (~isempty(sConeIndices))
        renderPatchArray(axesHandle, coneApertureShape, sConeApertureRadii, ...
            sConeAperturePositions, 3/4*0.9, edgeColor, lineWidth, lineStyle, faceAlpha);
    end

    % Show connections
    if (conesAreInputsToSingleRGC) && (~isempty(targetRGC))
        
        zLevels = [0.4 0.41];
        alpha = 0.4;
        cMap = brewermap(64, '*greys');
        contourLineColor = [0 0 0];
        cMosaic.semiTransparentContourPlot(axesHandle, xSupport, ySupport, ...
            rf, zLevels, cMap, alpha, contourLineColor);
        
        for k = 1:numel(coneIndices)
            theConePos = conePos(coneIndices(k),:);
            plot(axesHandle, [RGCpos(targetRGC,1) theConePos(1)], [RGCpos(targetRGC,2) theConePos(2)], 'k-', 'LineWidth', lineWidth*3);
            plot(axesHandle, [RGCpos(targetRGC,1) theConePos(1)], [RGCpos(targetRGC,2) theConePos(2)], '-', 'Color', [0.7 0.7 0.7], 'LineWidth', lineWidth*2);
        end
    end

    
    axis(axesHandle, 'equal')
    set(axesHandle, 'CLim', [0 1]);

end


function renderPatchArray(axesHandle, apertureShape, apertureRadii, rfCoords, ...
    faceColors, edgeColor, lineWidth, lineStyle, faceAlpha)

    conesNum = numel(apertureRadii);
    if (conesNum == 0)
        return;
    end
    
    verticesPerCone = numel(apertureShape.x);
    
    verticesList = zeros(verticesPerCone * conesNum, 2);
    facesList = [];
    
    if (numel(faceColors) == 1)
        colors = repmat(faceColors, [verticesPerCone*conesNum 1]);
    else
        colors = [];
    end
    
    for coneIndex = 1:conesNum
        idx = (coneIndex - 1) * verticesPerCone + (1:verticesPerCone);
        verticesList(idx, 1) = apertureShape.x*apertureRadii(coneIndex) + rfCoords(coneIndex,1);
        verticesList(idx, 2) = apertureShape.y*apertureRadii(coneIndex) + rfCoords(coneIndex,2);
        if ((numel(faceColors) == conesNum) && (conesNum > 1))
            colors = cat(1, colors, repmat(faceColors(coneIndex), [verticesPerCone 1]));
        end
        facesList = cat(1, facesList, idx);
    end

    S.Vertices = verticesList;
    S.Faces = facesList;
    S.FaceVertexCData = colors;
    S.FaceColor = 'flat';
    S.EdgeColor = edgeColor;
    S.FaceAlpha = faceAlpha;
    S.LineWidth = lineWidth;
    S.LineStyle = lineStyle;
    patch(S, 'Parent', axesHandle);
end
