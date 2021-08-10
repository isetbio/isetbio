% Method to visualize how the cone mosaic is tesselated by the RGC RF centers
function visualizeConeMosaicTesselation(obj, ...
    coneRFpositions, coneRFspacings, ...
    rgcRFpositions, rgcRFspacings, domain, ...
    varargin)

    p = inputParser;
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('showConnectedCones', true, @islogical);
    p.addParameter('fontSize', 16, @isscalar);
    p.addParameter('backgroundColor', [0.7 0.7 0.7]);
    p.addParameter('plotTitle', '', @ischar);
    p.addParameter('visualizationLimits', [], @(x)((isempty(x))||(numel(x)==4)));
    p.parse(varargin{:});
    
    figureHandle = p.Results.figureHandle;
    axesHandle = p.Results.axesHandle;
    visualizationLimits = p.Results.visualizationLimits;
    showConnectedCones = p.Results.showConnectedCones;
    fontSize = p.Results.fontSize;
    backgroundColor = p.Results.backgroundColor;
    plotTitle = p.Results.plotTitle;
    

    % Set figure size
    if (isempty(figureHandle))
        figureHandle = figure(); clf;
        set(figureHandle, 'Position', [10 10 700 700], 'Color', [1 1 1]);
        axesHandle = subplot('Position', [0.09 0.07 0.90 0.92]);
    else
        if (isempty(axesHandle))
            clf;
            set(figureHandle, 'Position', [10 10 700 700], 'Color', [1 1 1]);
            axesHandle = subplot('Position', [0.07 0.07 0.92 0.92]);
        end
        cla(axesHandle);
    end
    
    if (~isempty(visualizationLimits))
       idx = find( ...
           (rgcRFpositions(:,1) >= visualizationLimits(1)) & ...
           (rgcRFpositions(:,1) <= visualizationLimits(2)) & ...
           (rgcRFpositions(:,2) >= visualizationLimits(3)) & ...
           (rgcRFpositions(:,2) <= visualizationLimits(4)) );
       rgcRFpositions = rgcRFpositions(idx,:);
       rgcRFspacings = rgcRFspacings(idx);
       coneConnectivityMatrix = obj.coneConnectivityMatrix(:, idx);
    else
       coneConnectivityMatrix = obj.coneConnectivityMatrix;
    end
    
    if (isempty(rgcRFpositions))
        fprintf('No RGCs within visualizationLimits. Skipping visualization.\n');
        return;
    end
    
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
        
        
    % Spatial support
    nSamples = 200;
    xAxis = linspace(xRange(1), xRange(2), nSamples);
    yAxis = linspace(yRange(1), yRange(2), nSamples);
    [X,Y] = meshgrid(xAxis,yAxis);
    
    
    % Rendering params
    faceAlpha = 0.4;
    edgeAlpha = 0.8;
    zLevels = [0.3 1];
    whichLevelsToContour = [1];
        
    % Plot
    axis(axesHandle, 'equal');
    set(axesHandle, 'XLim', xRange, 'YLim', yRange, 'FontSize', 18);
    hold(axesHandle, 'on');
    
    % Plot them according to their ecc
    rgcRadialEcc = sum(rgcRFpositions.^2, 2);
    [~,sortedRGCindices] = sort(rgcRadialEcc);
    rgcsNum = numel(sortedRGCindices);
    
    for iRGC = 1:rgcsNum    
        
        RGCindex = sortedRGCindices(iRGC);
        
        % Find cones which are connected to this RGC
        connectivityVector = full(squeeze(coneConnectivityMatrix(:, RGCindex)));
        inputConeIDs = find(connectivityVector > 0.01);
        inputsNum = numel(inputConeIDs);
        
        if (inputsNum == 0)
            continue;
        end
        
        % Generate RF center outline based on cone inputs
        theRF = rfFromConnectivityMatrix(connectivityVector/max(connectivityVector), ...
            coneRFpositions, coneRFspacings, X,Y);
        
        C = contourc(xAxis, yAxis,theRF, zLevels);
        fillRFoutline(axesHandle, C, zLevels, whichLevelsToContour, faceAlpha, edgeAlpha);
    
        if (showConnectedCones)
            % Connected cones
            indicesOfConeInputsToThisRGC = find(connectivityVector>0);
            displayConnectedCones(axesHandle, indicesOfConeInputsToThisRGC, coneRFpositions);
        end
        
        if (mod(RGCindex-1,10) == 9)
            drawnow;
        end
        
    end % RGCindex
    
    % Display L-cones
    idx = find( ...
        (coneRFpositions(obj.inputConeMosaic.lConeIndices,1) >= xRange(1)) & ...
        (coneRFpositions(obj.inputConeMosaic.lConeIndices,1) <= xRange(2)) & ...
        (coneRFpositions(obj.inputConeMosaic.lConeIndices,2) >= yRange(1)) & ...
        (coneRFpositions(obj.inputConeMosaic.lConeIndices,2) <= yRange(2)));
    scatter(axesHandle,coneRFpositions(obj.inputConeMosaic.lConeIndices(idx),1), coneRFpositions(obj.inputConeMosaic.lConeIndices(idx),2), 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0.5 0.5]);
    
    % Display M-cones
    idx = find( ...
        (coneRFpositions(obj.inputConeMosaic.mConeIndices,1) >= xRange(1)) & ...
        (coneRFpositions(obj.inputConeMosaic.mConeIndices,1) <= xRange(2)) & ...
        (coneRFpositions(obj.inputConeMosaic.mConeIndices,2) >= yRange(1)) & ...
        (coneRFpositions(obj.inputConeMosaic.mConeIndices,2) <= yRange(2)));
    scatter(axesHandle,coneRFpositions(obj.inputConeMosaic.mConeIndices(idx),1), coneRFpositions(obj.inputConeMosaic.mConeIndices(idx),2), 'MarkerEdgeColor', [0 0.7 0], 'MarkerFaceColor', [0.5 0.9 0.5]);
    
    
    % Display S-cones
    idx = find( ...
        (coneRFpositions(obj.inputConeMosaic.sConeIndices,1) >= xRange(1)) & ...
        (coneRFpositions(obj.inputConeMosaic.sConeIndices,1) <= xRange(2)) & ...
        (coneRFpositions(obj.inputConeMosaic.sConeIndices,2) >= yRange(1)) & ...
        (coneRFpositions(obj.inputConeMosaic.sConeIndices,2) <= yRange(2)));
    scatter(axesHandle,coneRFpositions(obj.inputConeMosaic.sConeIndices(idx),1), coneRFpositions(obj.inputConeMosaic.sConeIndices(idx),2), 'MarkerEdgeColor', [0 0 1], 'MarkerFaceColor', [0.5 0.5 1.0]);
    
    
    % Finish plot
    xlabel(axesHandle, sprintf('space (%s)', domain));
    box(axesHandle, 'on');
    title(plotTitle);

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