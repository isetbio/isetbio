function visualize(obj, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('fillColor', [], @(x)(isempty(x)||((isvector(x))&&(numel(x) == 3))));
    p.parse(varargin{:});
    
    figureHandle = p.Results.figureHandle;
    axesHandle = p.Results.axesHandle;
    
    fillColor = p.Results.fillColor;
    if (isempty(fillColor))
        fillColor = [0.3 0.6 0.9];
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
    
    % Draw it
    hold(axesHandle, 'on');
    patch(axesHandle, roiOutline.x, roiOutline.y, ...
        fillColor, 'FaceAlpha', 0.5, 'EdgeColor', fillColor*0.5, 'LineWidth', 1.0);
    dx = 0.1*(max(roiOutline.x)-min(roiOutline.x));
    dy = 0.1*(max(roiOutline.y)-min(roiOutline.y));
    dx = max([dx dy]);
    axis(axesHandle, 'image');
    set(axesHandle, 'XLim', [min(roiOutline.x)-dx max(roiOutline.x)+dx], ...
                    'YLim', [min(roiOutline.y)-dx max(roiOutline.y)+dx]);
    xlabel(sprintf('%s', obj.geometryStruct.units));
    ylabel(sprintf('%s', obj.geometryStruct.units));
    set(axesHandle, 'FontSize', 14);
    grid(axesHandle, 'on');
    box(axesHandle, 'on');
end

