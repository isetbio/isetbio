function visualize(obj, varargin)
    p = inputParser;
    p.addParameter('domain', 'degrees', @(x)(ischar(x) && (ismember(x, {'degrees', 'microns'}))));
    p.addParameter('domainVisualizationLimits', [], @(x)((isempty(x))||(numel(x)==4)));
    p.addParameter('domainVisualizationTicks', [], @(x)(isempty(x)||(isstruct(x))));
    p.addParameter('covisualizeInputConeMosaic', false, @islogical);
    p.addParameter('covisualizeInputConeMosaicConnectivity', false, @islogical);
    p.addParameter('activation', []);
    p.addParameter('activationRange', [],@(x)((isempty(x))||(numel(x)==2)));
    p.addParameter('labelRetinalMeridians', false, @islogical);
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('fontSize', 16, @isscalar);
    p.addParameter('backgroundColor', [0.7 0.7 0.7]);
    p.addParameter('plotTitle', '', @ischar);
    p.parse(varargin{:});
    
    domain = p.Results.domain;
    domainVisualizationLimits = p.Results.domainVisualizationLimits;
    domainVisualizationTicks = p.Results.domainVisualizationTicks;
    covisualizeInputConeMosaic = p.Results.covisualizeInputConeMosaic;
    covisualizeInputConeMosaicConnectivity = p.Results.covisualizeInputConeMosaicConnectivity;
    activation = p.Results.activation;
    activationRange = p.Results.activationRange;
    labelRetinalMeridians = p.Results.labelRetinalMeridians;
    noXlabel = p.Results.noXLabel;
    noYlabel = p.Results.noYLabel;
    figureHandle = p.Results.figureHandle;
    axesHandle = p.Results.axesHandle;
    fontSize = p.Results.fontSize;
    backgroundColor = p.Results.backgroundColor;
    plotTitle = p.Results.plotTitle;
    
    % Determine displayed domain (degs or microns)
    switch (domain)
        case 'degrees'
            rfPositions = obj.rgcRFpositionsDegs;
            rfSpacings = obj.rgcRFspacingsDegs;
            coneRFpositions = obj.inputConeMosaic.coneRFpositionsDegs;
            rfProximityThreshold = 1/270;

        case 'microns'
            rfPositions = obj.coneRFpositionsMicrons;
            rfSpacings = obj.coneRFspacingsMicrons;
            coneRFPositions = obj.inputConeMosaic.coneRFpositionsMicrons;
            rfProximityThreshold = 1;

    end
    
    % Determine X,Y limits
    if (isempty(domainVisualizationLimits))
        xRange(1) = min(rfPositions(:,1));
        xRange(2) = max(rfPositions(:,1));
        if (xRange(2) == xRange(1))
            switch (domain)
                case 'degrees'
                    xRange = xRange(1) + 0.02*[-1 1];
                case 'microns'
                    xRange = xRange(1) + 2*[-1 1];
            end
        end
        yRange(1) = min(rfPositions(:,2));
        yRange(2) = max(rfPositions(:,2));
        if (yRange(2) == yRange(1))
            switch (domain)
                case 'degrees'
                    yRange = yRange(1) + 0.02*[-1 1];
                case 'microns'
                    yRange = yRange(1) + 2*[-1 1];
            end
        end
        xx = xRange(2)-xRange(1);
        yy = yRange(2)-yRange(1);
        xRange(1) = xRange(1)-xx*0.02;
        xRange(2) = xRange(2)+xx*0.02;
        yRange(1) = yRange(1)-yy*0.02;
        yRange(2) = yRange(2)+yy*0.02;
    else
        xRange(1) = domainVisualizationLimits(1);
        xRange(2) = domainVisualizationLimits(2);
        yRange(1) = domainVisualizationLimits(3);
        yRange(2) = domainVisualizationLimits(4);
    end
    
    % Set figure size
    if (isempty(figureHandle))
        figureHandle = figure(); clf;
        set(figureHandle, 'Position', [10 10 700 700], 'Color', [1 1 1]);
        axesHandle = subplot('Position', [0.05 0.05 0.94 0.94]);
    else
        if (isempty(axesHandle))
            figure(figureHandle);
            clf;
            set(figureHandle, 'Position', [10 10 700 700], 'Color', [1 1 1]);
            axesHandle = subplot('Position', [0.05 0.05 0.94 0.94]);
        end
        cla(axesHandle);
    end
    
    % Number of rgcs
    rgcsNum = numel(rfSpacings);
    
    
    % RF center shape (disk)
    if (rgcsNum > 10000)
        deltaAngle = 60;
    elseif (rgcsNum  > 5000)
        deltaAngle = 45;
    elseif (rgcsNum  > 1000)
       deltaAngle = 30;
    elseif (rgcsNum  > 500)
        deltaAngle = 20;
    elseif (rgcsNum  > 250)
        deltaAngle = 15;
    elseif (rgcsNum  > 100)
        deltaAngle = 10;
    else
        deltaAngle = 5;
    end
    iTheta = (0:deltaAngle:360) / 180 * pi;
    rfCenterShape.x = cos(iTheta);
    rfCenterShape.y = sin(iTheta);
    
    
    
    if (covisualizeInputConeMosaic)
        domainVisualizationTicks.x = [];
        domainVisualizationTicks.y = [];
        obj.inputConeMosaic.visualize(...
            'figureHandle', figureHandle, ...
            'axesHandle', axesHandle, ...
            'visualizedConeAperture', 'geometricArea', ...
            'visualizeConeApertureThetaSamples', 12, ...
            'domain', domain, ...
            'domainVisualizationLimits', domainVisualizationLimits, ...
            'domainVisualizationTicks', domainVisualizationTicks, ...
            'noXLabel', noXlabel, ...
            'noYlabel', noYlabel, ...
            'fontSize', fontSize, ...
            'backgroundColor', backgroundColor, ...
            'plotTitle', '');
    end
    
    hold(axesHandle, 'on');
    
    
    % Only visualize connectivity for RGCs within the domainVisualizationLimits
    indicesWithinVisualizationLimits = find(... 
            (rfPositions(:,1) >= xRange(1)) & ...
            (rfPositions(:,1) <= xRange(2)) & ...
            (rfPositions(:,2) >= yRange(1)) & ...
            (rfPositions(:,2) <= yRange(2)) ...
            );
        
    % Plot the RGC mosaic (either positions, or activation)
    if (isempty(activation))
        % Only visualize mRGC positions if we do not show connectivity
        %if (~covisualizeInputConeMosaicConnectivity)
            % Plot mRGC positions as white, semitransparent disks
            faceAlpha = 0.5;
            lineWidth = 1.0;
            edgeColor = [0.3 0.3 0.3];
            faceColor = [0.7 0.7 0.7];
            renderPatchArray(axesHandle, rfCenterShape, rfSpacings(indicesWithinVisualizationLimits)*0.5, ...
                rfPositions(indicesWithinVisualizationLimits,:), faceColor, edgeColor, lineWidth, faceAlpha);
       % end
    else
        % Plot mRGC activations
    end
    
    % Visualize connectivity
    if (covisualizeInputConeMosaicConnectivity)  
        for iRGC = 1:numel(indicesWithinVisualizationLimits)   
            RGCindex = indicesWithinVisualizationLimits(iRGC);
            % Find cones which are connected to this RGC
            connectivityVector = full(squeeze(obj.coneConnectivityMatrix(:, RGCindex)));
            inputConeIDs = find(connectivityVector > 1e-5);

            % Plot a line connecting this RGC to each of its input cones
            for iCone = 1:numel(inputConeIDs)
                plot(axesHandle, ...
                    [rfPositions(RGCindex,1) coneRFpositions(inputConeIDs(iCone), 1)], ...
                    [rfPositions(RGCindex,2) coneRFpositions(inputConeIDs(iCone), 2)], ...
                    'k-', 'MarkerFaceColor', [0 0 0], 'MarkerSize', 5, 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.0);
            end
            % Plot a dot on each of the input cones
            plot(axesHandle, ...
                    coneRFpositions(inputConeIDs, 1), ...
                    coneRFpositions(inputConeIDs, 2), ...
                    'ko', 'MarkerSize', 6,'LineWidth', 1.0);
        end
    end
        
    hold(axesHandle, 'off');
    
    % Finalize plot
    set(axesHandle, 'Color', backgroundColor);
    axis(axesHandle, 'xy');
    axis(axesHandle, 'equal');
    set(axesHandle, 'XLim', xRange, 'YLim', yRange, 'CLim', [0 1], 'FontSize', fontSize);        
   
    if (isempty(domainVisualizationTicks))
        xo = (xRange(1)+xRange(2))/2;
        xx = xRange(2)-xRange(1);
        yo = (yRange(1)+yRange(2))/2;
        yy = yRange(2)-yRange(1);
        ticksX = xo + xx*0.5*[-0.75 0 0.75];
        ticksY = yo + yy*0.5*[-0.75 0 0.75];
        
        if (xx > 10)
            domainVisualizationTicks.x = round(ticksX);
        elseif (xx > 1)
            domainVisualizationTicks.x = round(ticksX*100)/100;
        else
            domainVisualizationTicks.x = round(ticksX*1000)/1000;
        end
        if (yy > 10)
            domainVisualizationTicks.y = round(ticksY);
        elseif (yy > 1)
            domainVisualizationTicks.y = round(ticksY*100)/100;
        else
            domainVisualizationTicks.y = round(ticksY*1000)/1000;
        end

    end
    
    set(axesHandle, 'XTick', domainVisualizationTicks.x, ...
                    'YTick', domainVisualizationTicks.y);
     
    box(axesHandle, 'on');
    set(figureHandle, 'Color', [1 1 1]);
    
    switch (domain)
        case 'degrees'
            if (~noXlabel)
                if (labelRetinalMeridians)
                    if (strcmp(obj.whichEye, 'right eye'))
                        leftMeridianName = 'temporal retina';
                        rightMeridianName = 'nasal retina';
                    else
                        leftMeridianName = 'nasal retina';
                        rightMeridianName = 'temporal retina';
                    end
                    xlabel(axesHandle, sprintf('\\color{red}%s    \\color{black} retinal space (degrees)    \\color[rgb]{0 0.7 0} %s', ...
                            leftMeridianName, rightMeridianName));
                else
                    xlabel(axesHandle, 'space (degrees)');
                end
            end
            if (~noYlabel)
                if (labelRetinalMeridians)
                    ylabel(axesHandle, sprintf('%s  < = = = = = |     space (degrees)    | = = = = =  > %s', ...
                        'superior retina', 'inferior retina'));
                else
                    ylabel(axesHandle, 'space (degrees)');
                end
            end
            set(axesHandle, 'XTickLabel', sprintf('%1.1f\n', domainVisualizationTicks.x), ...
                            'YTickLabel', sprintf('%1.1f\n', domainVisualizationTicks.y));
        case 'microns'
            if (~noXlabel)
                if (labelRetinalMeridians)
                    if (strcmp(obj.whichEye, 'right eye'))
                        leftMeridianName = '(temporal)';
                        rightMeridianName = '(nasal)';
                    else
                        leftMeridianName = '(nasal)';
                        rightMeridianName = '(temporal)';
                    end
                    xlabel(axesHandle, sprintf('\\color{red}%s    \\color{black} retinal space (microns)    \\color[rgb]{0 0.7 0} %s', ...
                            leftMeridianName, rightMeridianName));
                else
                    xlabel(axesHandle, 'retinal space (microns)');
                end
            end
            if (~noYlabel)
                if (labelRetinalMeridians)
                    upperMeridianName = '(inferior)';
                    lowerMeridianName = '(superior)';
                    ylabel(axesHandle, sprintf('\\color{blue}%s    \\color{black} retinal space (microns)    \\color[rgb]{0.6 0.6 0.4} %s', ...
                            lowerMeridianName, upperMeridianName));
                else
                    ylabel(axesHandle, 'space (microns)');
                end
            end
            set(axesHandle, 'XTickLabel', sprintf('%d\n', domainVisualizationTicks.x), ...
                            'YTickLabel', sprintf('%d\n', domainVisualizationTicks.y));
    end
    
    if (~isempty(plotTitle))
        title(axesHandle,plotTitle);
    end
    
end

function renderPatchArray(axesHandle, rfShape, rfRadii, rfCoords, ...
    faceColor, edgeColor, lineWidth, faceAlpha)

    rfsNum = numel(rfRadii);
    if (rfsNum == 0)
        return;
    end
    
    verticesPerRF = numel(rfShape.x);
    verticesList = zeros(verticesPerRF * rfsNum, 2);
    facesList = [];
    
   % if (numel(faceColors) == 1)
   %     colors = repmat(faceColors, [verticesPerRF*rfsNum 1]);
   % else
   %     colors = [];
   % end
    
    for rfIndex = 1:rfsNum
        idx = (rfIndex - 1) * verticesPerRF + (1:verticesPerRF);
        verticesList(idx, 1) = rfShape.x*rfRadii(rfIndex) + rfCoords(rfIndex,1);
        verticesList(idx, 2) = rfShape.y*rfRadii(rfIndex) + rfCoords(rfIndex,2);
%         if (numel(faceColors) == rfsNum)
%             colors = cat(1, colors, repmat(faceColors(rfIndex), [verticesPerRF 1]));
%         end
        facesList = cat(1, facesList, idx);
    end

    S.Vertices = verticesList;
    S.Faces = facesList;
   % S.FaceVertexCData = colors;
    S.FaceColor = faceColor;
    S.EdgeColor = edgeColor;
    S.FaceAlpha = faceAlpha;
    S.LineWidth = lineWidth;
    patch(S, 'Parent', axesHandle);
end
