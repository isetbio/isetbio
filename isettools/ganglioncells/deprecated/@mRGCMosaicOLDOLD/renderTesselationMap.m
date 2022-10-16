function renderTesselationMap(figNo, axesHandle, coneRFpositions, coneRFspacings, coneTypes, rgcRFpositions, rgcRFspacings, coneConnectivityMatrix, domain)

    % xRange
    [m, idx] = min(rgcRFpositions(:,1));
    xRange(1) = m - rgcRFspacings(idx);
    [m, idx] = max(rgcRFpositions(:,1));
    xRange(2) = m + rgcRFspacings(idx);

    % yRange
    [m, idx] = min(rgcRFpositions(:,2));
    yRange(1) = m - rgcRFspacings(idx);
    [m, idx] = max(rgcRFpositions(:,2));
    yRange(2) = m + rgcRFspacings(idx);
            
    
    
    % Reset stats
    stats.rgcsNumWithNonZeroInputs = 0;
    stats.coneInputsPerRGC = zeros(1,100);
    stats.mixedInputRGCs = 0;
    
    % Spatial support
    nSamples = 200;
    xAxis = linspace(xRange(1), xRange(2), nSamples);
    yAxis = linspace(yRange(1), yRange(2), nSamples);
    [X,Y] = meshgrid(xAxis,yAxis);
    
    faceAlpha = 0.4;
    edgeAlpha = 0.8;
    zLevels = [0.3 1];
    whichLevelsToContour = [1];
    showConnectedCones = true;
        
    % Plot
    if (isempty(figNo))
        hFig = figure(); clf;
        axesHandle = subplot('Position', [0.05 0.06 0.93 0.93]);
        set(hFig, 'Position', [10 10 1100 1100], 'Color', [1 1 1]);
    else
        hFig = figure(figNo); 
        if (isempty(axesHandle))
            clf;
            set(hFig, 'Position', [10 10 1100 1100], 'Color', [1 1 1]);
            axesHandle = subplot('Position', [0.05 0.06 0.93 0.93]);
        else
            set(hFig, 'Position', [10 10 2100 1100], 'Color', [1 1 1]);
            if (axesHandle == 1)
                axesHandle = subplot('Position', [0.05 0.06 0.45 0.93]);
            else
                axesHandle = subplot('Position', [0.52 0.06 0.45 0.93]);
            end
        end
    end

    
    
    hold(axesHandle, 'on');
    
    rgcsNum = size(rgcRFpositions ,1);
    for RGCindex = 1:rgcsNum    
        
        % Find number of cone inputs
        connectivityVector = full(squeeze(coneConnectivityMatrix(:, RGCindex)));
        inputIDs = find(connectivityVector > 0.01);
        inputsNum = numel(inputIDs);
        if (inputsNum == 0)
            fprintf(2, 'RGC %d has zero inputs!!!\n', RGCindex);
            continue;
        end
        
        % Update stats
        stats.rgcsNumWithNonZeroInputs = stats.rgcsNumWithNonZeroInputs + 1;
        stats.coneInputsPerRGC(inputsNum) = stats.coneInputsPerRGC(inputsNum) + 1;
        inputTypesNum = numel(unique(coneTypes(inputIDs)));
        if (inputTypesNum>1)
            stats.mixedInputRGCs = stats.mixedInputRGCs + 1;
        end
        
        % Generate RF center outline based on cone inputs
        theRF = rfFromConnectivityMatrix(connectivityVector/max(connectivityVector), ...
            coneRFpositions, coneRFspacings, X,Y);
        
        C = contourc(xAxis, yAxis,theRF, zLevels);
        fillRFoutline(axesHandle, C, zLevels, whichLevelsToContour, faceAlpha, edgeAlpha);
    
        indicesOfConeInputsToThisRGC = find(connectivityVector>0);
        
        if (showConnectedCones)
            % Connected cones
            displayConnectedCones(axesHandle, indicesOfConeInputsToThisRGC, coneRFpositions);
        end
        
        if (mod(RGCindex-1,10) == 9)
            drawnow;
        end
        
    end % RGCindex
    
    % Display cones
    LconeIndices = find(coneTypes == mRGCmosaic.LCONE_ID);
    MconeIndices = find(coneTypes == mRGCmosaic.MCONE_ID);
    SconeIndices = find(coneTypes == mRGCmosaic.SCONE_ID);
    scatter(axesHandle,coneRFpositions(LconeIndices,1), coneRFpositions(LconeIndices,2), 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0.5 0.5]);
    scatter(axesHandle,coneRFpositions(MconeIndices,1), coneRFpositions(MconeIndices,2), 'MarkerEdgeColor', [0 0.7 0], 'MarkerFaceColor', [0.5 0.9 0.5]);
    scatter(axesHandle,coneRFpositions(SconeIndices,1), coneRFpositions(SconeIndices,2), 'MarkerEdgeColor', [0 0 1], 'MarkerFaceColor', [0.5 0.5 1.0]);
    
    % Finish plot
    axis(axesHandle, 'equal');
    set(axesHandle, 'XLim', xRange, 'YLim', yRange, 'FontSize', 18);
    xlabel(axesHandle, sprintf('space (%s)', domain));
    box(axesHandle, 'on');
    set(hFig, 'Name', sprintf('cones:%d, mRGCs:%d \nSingle cone RF: %2.1f%%, 2-cone:  %2.1f%%, 3-cone:  %2.1f%%, 4-cone:  %2.1f%%, 5-cone:  %2.1f%%, 6+ cones:  %2.1f%%', ...
        size(coneConnectivityMatrix,1), size(coneConnectivityMatrix,2), ...
        100*stats.coneInputsPerRGC(1)/sum(stats.coneInputsPerRGC), ...
        100*stats.coneInputsPerRGC(2)/sum(stats.coneInputsPerRGC), ...
        100*stats.coneInputsPerRGC(3)/sum(stats.coneInputsPerRGC), ...
        100*stats.coneInputsPerRGC(4)/sum(stats.coneInputsPerRGC), ...
        100*stats.coneInputsPerRGC(5)/sum(stats.coneInputsPerRGC), ...
        100*sum(stats.coneInputsPerRGC(6:end))/sum(stats.coneInputsPerRGC) ...
        ));
end



function [theRF, theRFpos] = rfFromConnectivityMatrix(connectivityVector, conePositions, coneSpacings,  X,Y)
    
    theRF = [];
    connectedConeIDs = find(connectivityVector>0);
    flatTopZ = 0.45;

    theRFpos = mean(conePositions(connectedConeIDs,:),1);
    
    for k = 1:numel(connectedConeIDs)
        coneIndex = connectedConeIDs(k);
        cP = squeeze(conePositions(coneIndex,:));
        coneSigma = coneSpacings(coneIndex)/3;
        coneProfile = exp(-0.5*((X-cP(1))/coneSigma).^2) .* exp(-0.5*((Y-cP(2))/coneSigma).^2);
        coneProfile(coneProfile>=flatTopZ) = flatTopZ;
        if (isempty(theRF))
            theRF = coneProfile * connectivityVector(coneIndex);
        else
            theRF = theRF + coneProfile * connectivityVector(coneIndex);
        end 
    end
    if (~isempty(theRF))
        theRF = theRF/max(theRF(:));
    end
end

function  fillRFoutline(theAxes, C, zLevels, whichLevelsToContour, faceAlpha, edgeAlpha)
    k = 1;
    contoursNum = 0;
    while k < size(C,2)
        level = C(1,k);
        points = C(2,k);
        if  (level ~= zLevels(whichLevelsToContour(1))) && (~ismember(level, zLevels(whichLevelsToContour)))
            % skip this contour
            k = k+points+1;
            continue;
        end
        
        xOutline = C(1,k+(1:points));
        yOutline = C(2,k+(1:points));
        
        faceColor = [0.6 0.6 0.5]-level*0.05;
        edgeColor = [0.2 0.2 0.2];
        patchContour(theAxes, xOutline, yOutline, faceColor, edgeColor, faceAlpha, edgeAlpha);

        k = k+points+1;
        contoursNum = contoursNum + 1;
    end
end

function patchContour(theAxes, xRGCEnsembleOutline, yRGCEnsembleOutline, faceColor, edgeColor, faceAlpha, edgeAlpha)
    v = [xRGCEnsembleOutline(:) yRGCEnsembleOutline(:)];
    f = 1:numel(xRGCEnsembleOutline);
    patch(theAxes, 'Faces', f, 'Vertices', v, 'FaceColor', faceColor, ...
            'FaceAlpha', faceAlpha, 'EdgeColor', edgeColor, ... 
           'EdgeAlpha', edgeAlpha, 'LineWidth', 1.5);
end

function displayConnectedCones(ax, indicesOfConeInputs, conePositionsMicrons)
    % Polygon connecting input cones
    xx = conePositionsMicrons(indicesOfConeInputs,1);
    yy = conePositionsMicrons(indicesOfConeInputs,2);
    xo = mean(xx);
    yo = mean(yy);
    dx = xx-xo;
    dy = yy-yo;
    [~,idx] = sort(unwrap(atan2(dy,dx)));
    xx = xx(idx);
    yy = yy(idx);
    xx(end+1) = xx(1);
    yy(end+1) = yy(1);
    plot(ax, xx,yy, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
end