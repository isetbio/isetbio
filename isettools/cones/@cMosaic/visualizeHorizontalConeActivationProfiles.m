function figureHandle = visualizeHorizontalConeActivationProfiles(obj, theConeMosaicResponse, coneTypesToVisualize, ...
    horizontalSliceYcoordDegs, horizontalSliceThicknessDegs, maxResponse, varargin)
% Superimpose cone activation profiles along a horizontal slice on top of the 2D cMosaic activation pattern
%
% Syntax:
%   hFig = visualizeHorizontalConeActivationProfiles(obj, ...
%    theConeMosaicResponse, coneTypesToVisualize, ...
%    horizontalSliceYcoordDegs, horizontalSliceThicknessDegs, ...
%    maxResponse, varargin);
%
% Description:
%    Superimpose cone activation profiles along a horizontal slice on top
%    of the 2D cMosaic activation pattern
%
% Inputs:
%    obj                            - A @cMosaic object
%    theConeMosaicResponse          - A @cMosaic activation pattern
%    coneTypesToVisualize           - An array of cone types to visualize,
%                                     e.g., [cMosaic.LCONE_ID, cMosaic.MCONE_ID, cMosaic.SCONE_ID]
%    horizontalSliceYcoordDegs      - Y-coord of slice
%    horizontalSliceThicknessDegs   - Y-thickness of slice
%    maxResponse                    - The maximum value of the response
%
% Outputs:
%    figureHandle                   - Handle of generated figure
%
% Optional key/value pairs:
%    'figureHandle'                 - Figure handle on which to render
%    'axesHandle'                   - Axes handle on which to render
%    'onTopOfConeMosaicActivation'  - Render plot on top of the coneMosaicActivation pattern
%    'activationProfileScalingDegs' - Y-axis scaling (in degrees) of
%                                     activation profiles (this is elevant only if
%                                     only if 'onTopOfConeMosaicActivation' is true)

    % Parse input
    p = inputParser;
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('clearAxesBeforeDrawing', false, @islogical);
    p.addParameter('onTopOfConeMosaicActivation', true, @islogical);
    p.addParameter('activationProfileScalingDegs', 0.5, @isnumeric);
    p.addParameter('fontSize', 16, @isscalar);
    p.parse(varargin{:});

    figureHandle = p.Results.figureHandle;
    axesHandle = p.Results.axesHandle;
    renderOnTopOfConeMosaicActivation = p.Results.onTopOfConeMosaicActivation;
    activationProfileScalingDegs = p.Results.activationProfileScalingDegs;
    clearAxesBeforeDrawing = p.Results.clearAxesBeforeDrawing;
    fontSize = p.Results.fontSize;

    % Visualized L-cone indices and x-coords
    if (ismember(cMosaic.LCONE_ID, coneTypesToVisualize))
        [theVisualizedLConeIndices, theROI] = ...
            obj.coneIndicesOfCertainTypeWithinHorizontalSlice(horizontalSliceYcoordDegs, horizontalSliceThicknessDegs, cMosaic.LCONE_ID);
    else
        theVisualizedLConeIndices = [];
    end

    % Visualized M-cone indices and x-coords
    if (ismember(cMosaic.MCONE_ID, coneTypesToVisualize))
        theVisualizedMConeIndices = ...
            obj.coneIndicesOfCertainTypeWithinHorizontalSlice(horizontalSliceYcoordDegs, horizontalSliceThicknessDegs, cMosaic.MCONE_ID);
    else
        theVisualizedMConeIndices = [];
    end


    % Visualized S-cone indices and x-coords
    if (ismember(cMosaic.SCONE_ID, coneTypesToVisualize))
        theVisualizedSConeIndices = ...
            obj.coneIndicesOfCertainTypeWithinHorizontalSlice(horizontalSliceYcoordDegs, horizontalSliceThicknessDegs, cMosaic.SCONE_ID);
    else
        theVisualizedSConeIndices = [];
    end

    visualizationLimits = [...
        obj.eccentricityDegs(1) - obj.sizeDegs(1)*0.51 ...
        obj.eccentricityDegs(1) + obj.sizeDegs(1)*0.51 ...
        obj.eccentricityDegs(2) - obj.sizeDegs(2)*0.51 ...
        obj.eccentricityDegs(2) + obj.sizeDegs(2)*0.51 ...
        ];
    visualizationTicks = struct('x', -20:0.5:20, 'y', -20:0.5:20);

    % Set figure size
    if (isempty(figureHandle))
        figureHandle = figure();clf;
        set(figureHandle, 'Position', [10 10 1900 1050], 'Color', [0 0 0]);
        if (renderOnTopOfConeMosaicActivation)
            axesHandle = subplot('Position', [0.05 0.05 0.94 0.94]);
        else
            axesHandle = subplot('Position', [0.06 0.09 0.93 0.9]);
        end
    else
        if (isempty(axesHandle))
            figure(figureHandle);
            clf;
            set(figureHandle, 'Position', [10 10 1900 1050], 'Color', [0 0 0]);
            if (renderOnTopOfConeMosaicActivation)
                axesHandle = subplot('Position', [0.05 0.05 0.94 0.94]);
            else
                axesHandle = subplot('Position', [0.06 0.09 0.93 0.9]);
            end
        end

        if (clearAxesBeforeDrawing)
            cla(axesHandle);
        end
    end

    if (renderOnTopOfConeMosaicActivation)
        % Plot the cone mosaic response
        obj.visualize(...
            'figureHandle', figureHandle, ...
            'axesHandle', axesHandle, ...
            'activation', theConeMosaicResponse, ...
            'verticalActivationColorBarInside', true, ...
            'backgroundColor', [0.2 0.2 0.2], ...
            'domainVisualizationLimits', visualizationLimits, ...
            'domainVisualizationTicks', visualizationTicks);
        hold(axesHandle, 'on');

        % ROI over which we visualize cone responses
        theROI.visualize(...
            'figureHandle', figureHandle, ...
            'axesHandle', axesHandle, ...
            'xLims', visualizationLimits(1:2), ...
            'yLims', visualizationLimits(3:4), ...
            'fillColor', [0.1 1.0 0.1 0.6]);
    

        % Y-offset
        yOffset = horizontalSliceThicknessDegs*0.5+0.1;
    
        % Render visualization background
        renderVisualizationBackground(axesHandle, obj.eccentricityDegs, obj.sizeDegs, ...
            horizontalSliceYcoordDegs, yOffset, activationProfileScalingDegs);
    

        % Render L-cone responses
        if (~isempty(theVisualizedLConeIndices))
            superimposeConeResponses(axesHandle, squeeze(obj.coneRFpositionsDegs(theVisualizedLConeIndices,1)), theVisualizedLConeIndices, ...
                theConeMosaicResponse, horizontalSliceYcoordDegs+yOffset, activationProfileScalingDegs, maxResponse, [1.0 0.2 0.4]);
        end
    
        % Render M-cone responses
        if (~isempty(theVisualizedMConeIndices))
            superimposeConeResponses(axesHandle, squeeze(obj.coneRFpositionsDegs(theVisualizedMConeIndices,1)), theVisualizedMConeIndices, ...
                theConeMosaicResponse, horizontalSliceYcoordDegs+yOffset, activationProfileScalingDegs, maxResponse, [0.1 1.0 0.6]);
        end
    
    
        % Render S-cone responses
        if (~isempty(theVisualizedSConeIndices))
            superimposeConeResponses(axesHandle, squeeze(obj.coneRFpositionsDegs(theVisualizedSConeIndices,1)), theVisualizedSConeIndices, ...
                theConeMosaicResponse, horizontalSliceYcoordDegs+yOffset, activationProfileScalingDegs, maxResponse, [0.6 0.1 1.0]);
        end
    
        % Finalize plot
        hold(axesHandle, 'off')
        set(figureHandle, 'Color', [0 0 0]);
        set(axesHandle, 'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], 'LineWidth', 1.5, 'FontSize', fontSize);
        grid(axesHandle, 'off')

    else
        % Render L-cone responses
        if (~isempty(theVisualizedLConeIndices))
            renderConeResponses(axesHandle, squeeze(obj.coneRFpositionsDegs(theVisualizedLConeIndices,1)), theVisualizedLConeIndices, ...
                theConeMosaicResponse,  maxResponse, [1.0 0.2 0.4]);
            hold(axesHandle, 'on')
        end
    
        % Render M-cone responses
        if (~isempty(theVisualizedMConeIndices))
            renderConeResponses(axesHandle, squeeze(obj.coneRFpositionsDegs(theVisualizedMConeIndices,1)), theVisualizedMConeIndices, ...
                theConeMosaicResponse, maxResponse, [0.1 1.0 0.6]);
            hold(axesHandle, 'on')
        end
    
    
        % Render S-cone responses
        if (~isempty(theVisualizedSConeIndices))
            renderConeResponses(axesHandle, squeeze(obj.coneRFpositionsDegs(theVisualizedSConeIndices,1)), theVisualizedSConeIndices, ...
                theConeMosaicResponse, maxResponse, [0.6 0.1 1.0]);
        end

        hold(axesHandle, 'off');

        set(figureHandle, 'Color', [1 1 1]);
        set(axesHandle, 'LineWidth', 1.5, 'FontSize', fontSize);
        grid(axesHandle, 'on');
        xlabel(axesHandle, 'degs');
        ylabel(axesHandle, sprintf('excitations/%2.0fmsec', obj.integrationTime*1000));
    end


end


function renderVisualizationBackground(ax, mosaicEccDegs, mosaicSizeDegs, horizontalSliceYcoordDegs, yOffsetDegs, visualizedResponseScalingDegs)
    xx = mosaicEccDegs(1)+0.5*mosaicSizeDegs(1)*[-1 1 1 -1 -1];
    ymin = -visualizedResponseScalingDegs*0.1;
    ymax =  visualizedResponseScalingDegs*1.1;
    yy = horizontalSliceYcoordDegs + yOffsetDegs + [ymin ymin ymax ymax ymin];
    yGridMid = horizontalSliceYcoordDegs + yOffsetDegs+ visualizedResponseScalingDegs*0.5;
    vertices = [xx(:) yy(:)];
    patch(ax, 'Faces',1:size(vertices,1),'Vertices',vertices, 'EdgeColor', [1 1 0], 'EdgeAlpha', 0.7, 'FaceColor', [0.1 0.1 0.1], 'FaceAlpha', 0.85);
    plot(ax, xx, xx*0+yGridMid, 'y:','LineWidth', 0.5);
end


function superimposeConeResponses(ax, theVisualizedConeXcoords, theVisualizedConeIndices, ...
    theConeMosaicResponse, horizontalSliceYcoord, visualizedResponseScalingDegs, maxResponse, coneColor)

    theVisualizedConeResponses = squeeze(theConeMosaicResponse(1,1,theVisualizedConeIndices));
    scatter(ax, theVisualizedConeXcoords, horizontalSliceYcoord + theVisualizedConeResponses/maxResponse*visualizedResponseScalingDegs, 100, ...
        'MarkerEdgeColor', coneColor, 'MarkerFaceColor', coneColor, 'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 1.0, 'LineWidth', 1.5);
end

function renderConeResponses(ax, theVisualizedConeXcoords, theVisualizedConeIndices, ...
    theConeMosaicResponse, maxResponse, coneColor)

    theVisualizedConeResponses = squeeze(theConeMosaicResponse(1,1,theVisualizedConeIndices));
    scatter(ax, theVisualizedConeXcoords, theVisualizedConeResponses, 100, ...
        'MarkerEdgeColor', coneColor, 'MarkerFaceColor', coneColor, 'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 1.0, 'LineWidth', 1.5);
    set(ax, 'YLim', [0 maxResponse])
end
