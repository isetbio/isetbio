function visualizeActivationMap(obj, axesHandle, activation, varargin)

    % Parse the input
    p = inputParser;
    p.addRequired('axesHandle', @ishandle);
    p.addRequired('activation', @isnumeric);
    p.addParameter('domain', 'degs', @(x)(ischar(x) && (ismember(x, {'degs', 'microns'}))));
    p.addParameter('signalRange', [], @isnumeric);
    p.addParameter('visualizedFOV', [], @(x)(isnumeric(x) && (numel(x) == 2)));
    p.addParameter('backgroundColor', [0 0 0], @(x)(isnumeric(x) && (numel(x) == 3)));
    
    p.addParameter('fontSize', 18, @isnumeric);
    p.parse(axesHandle, activation, varargin{:});
    
    domain = p.Results.domain;
    signalRange = p.Results.signalRange;
    visualizedFOV = p.Results.visualizedFOV;
    backgroundColor = p.Results.backgroundColor;
    fontSize = p.Results.fontSize;
    
    if (isempty(signalRange))
        signalRange = [min(activation(:)) max(activation(:))];
    end
    
    switch (domain)
        case 'degs'
            rgcRFpositions = obj.rgcRFpositionsDegs;
            rgcRFspacings = obj.rgcRFspacingsDegs;
        case 'microns'
            rgcRFpositions = obj.rgcRFpositionsMicrons;
            rgcRFspacings = obj.rgcRFspacingsMicrons;
    end

    
    if (isempty(visualizedFOV))
        xRange = [min(rgcRFpositions(:,1)) max(rgcRFpositions(:,1))];
        yRange = [min(rgcRFpositions(:,2)) max(rgcRFpositions(:,2))];
    else
        xRange = mean(rgcRFpositions(:,1)) + visualizedFOV(1)*[-0.5 0.5];
        yRange = mean(rgcRFpositions(:,2)) + visualizedFOV(2)*[-0.5 0.5];
    end
    
    rgcsNum = size(rgcRFpositions,1);
    verticesNum = 10;
    deltaAngle = 360/(verticesNum-1);
    angles = (0:deltaAngle:360)';
    xOutline = cosd(angles);
    yOutline = sind(angles);
    verticesList = zeros(rgcsNum*verticesNum, 2);
    facesList = [];
    colors = [];
    
    for iRGC = 1: rgcsNum
        idx = (iRGC - 1) * verticesNum + (1:verticesNum);    
        faceColor = (activation(iRGC)-signalRange(1))/(signalRange(2)-signalRange(1)) * [1 1 1];
        radius = rgcRFspacings(iRGC) * 0.4;
        verticesList(idx, 1) = rgcRFpositions(iRGC,1) + xOutline * radius;
        verticesList(idx, 2) = rgcRFpositions(iRGC,2) + yOutline * radius;
        facesList = cat(1, facesList, idx);
        colors = cat(1, colors, repmat(faceColor, [verticesNum 1]));
    end

    S.Vertices = verticesList;
    S.Faces = facesList;
    S.FaceVertexCData = colors;
    S.FaceColor = 'flat';
    S.EdgeColor = [0 0 0];
    S.LineWidth = 1.0;
    patch(S, 'Parent', axesHandle);
    
    % Finalize plot
    axis(axesHandle, 'equal');
    set(axesHandle, 'XLim', xRange, 'YLim', yRange, 'XTick',-10:0.25:10, 'YTick', -10:0.25:10, 'Color', backgroundColor, 'FontSize', fontSize);
end
