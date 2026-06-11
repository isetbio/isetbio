function visualize(obj, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('xLims', [], @(x)(isempty(x)||((isvector(x))&&(numel(x) == 2))));
    p.addParameter('yLims', [], @(x)(isempty(x)||((isvector(x))&&(numel(x) == 2))));
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.addParameter('fillColor', [], @(x)(isempty(x)||((isvector(x))&&((numel(x) == 3)||(numel(x)==4)))));
    p.addParameter('fontSize', 14, @isscalar);
    p.parse(varargin{:});
    
    figureHandle = p.Results.figureHandle;
    axesHandle = p.Results.axesHandle;
    fillColor = p.Results.fillColor;
    xLims = p.Results.xLims;
    yLims = p.Results.yLims;
    fontSize = p.Results.fontSize;
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;


    if (isempty(fillColor))
        fillColor = [0.3 0.6 0.9];
        faceAlpha = 0.5;
    else
        if (numel(fillColor) == 4)
            faceAlpha = fillColor(4);
            fillColor = fillColor(1:3);
        else
            faceAlpha = 0.5;
        end
    end
    
    

    % Set figure size
    if (isempty(figureHandle))
        figureHandle = figure(); clf;
        set(figureHandle, 'Position', [10 10 700 700], 'Color', [1 1 1]);
        axesHandle = subplot('Position', [0.09 0.07 0.85 0.90]);
    else
        if (isempty(axesHandle))
            figure(figureHandle);
            clf;
            set(figureHandle, 'Position', [10 10 700 700], 'Color', [1 1 1]);
            axesHandle = subplot('Position', [0.09 0.07 0.85 0.90]);
        end
        %cla(axesHandle);
    end
    
    % Generate the outline
    roiOutline = obj.outline();
    
    dx = 0.1*(max(roiOutline.x)-min(roiOutline.x));
    dy = 0.1*(max(roiOutline.y)-min(roiOutline.y));
    dx = max([dx dy]);
    if (isempty(xLims))
        xLims = [min(roiOutline.x)-dx max(roiOutline.x)+dx];
    end

    if (isempty(yLims))
        yLims = [min(roiOutline.y)-dx max(roiOutline.y)+dx];
    end

    % Draw it
    hold(axesHandle, 'on');
    patch(axesHandle, roiOutline.x, roiOutline.y, ...
        fillColor, 'FaceAlpha', faceAlpha, 'EdgeColor', fillColor*0.5, 'LineWidth', 1.0);
    
    
    axis(axesHandle, 'image');
    set(axesHandle, 'XLim', xLims, ...
                    'YLim', yLims);
    if (~noXLabel)
        xlabel(sprintf('%s', obj.geometryStruct.units));
    end
    if (~noYLabel)
        ylabel(sprintf('%s', obj.geometryStruct.units));
    end
    set(axesHandle, 'FontSize', fontSize);
    grid(axesHandle, 'on');
    box(axesHandle, 'on');
end

