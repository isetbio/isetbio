function renderMultiSTF(ax, sfSupportCPD, theSTFs, theShadedSTF, theSTFphases, theSTFColors, theLegends, plotTitle, ff, varargin)
    p = inputParser;
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.addParameter('noXTickLabel', false, @islogical);
    p.addParameter('noYTickLabel', false, @islogical);
    p.addParameter('visualizedSpatialFrequencyRange', [], @(x)(isempty(x)||(numel(x)==2)));
    p.parse(varargin{:});
    
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;
    noXTickLabel = p.Results.noXTickLabel;
    noYTickLabel = p.Results.noYTickLabel;
    visualizedSpatialFrequencyRange = p.Results.visualizedSpatialFrequencyRange;

    yyaxis(ax, 'right')
    if (~isempty(theShadedSTF))
        y1 = theShadedSTF;
        y2 = 0*y1;
        faceColor = 0.5*(squeeze(theSTFColors(1,:))) + [0.5 0.5 0.5];
        edgeColor = squeeze(theSTFColors(1,:));
        faceAlpha = 0.3;
        lineWidth = 1.0;
        lineStyle = '-';
        shadedAreaBetweenTwoLines(ax, sfSupportCPD, y1,y2, ...
            faceColor, edgeColor, faceAlpha, lineWidth, lineStyle);
        hold(ax, 'on');
    end

    for iSTF = 1:size(theSTFs,2)
        pH(iSTF) = plot(ax, sfSupportCPD, squeeze(theSTFs(:,iSTF)), 'o-', 'LineWidth', ff.lineWidth*2, ...
            'MarkerEdgeColor', squeeze(theSTFColors(iSTF,:)), ...
            'MarkerFaceColor', 0.3*squeeze(theSTFColors(iSTF,:)) + 0.7*[1 1 1], ...
            'Color', squeeze(theSTFColors(iSTF,:)), ...
            'MarkerSize', ff.markerSize, ...
            'LineWidth', ff.lineWidth);
        hold(ax, 'on');
    end

    axis(ax, 'square');
    
    % ticks and grids
    if (isempty(visualizedSpatialFrequencyRange))
        visualizedSpatialFrequencyRange(1) = min(sfSupportCPD);
        visualizedSpatialFrequencyRange(2) = max(sfSupportCPD);
    end
    if (visualizedSpatialFrequencyRange(2) <= 40)
        xTicks = [0.01 0.03 0.1 0.3 1 3 10 20 30 40];
    elseif (visualizedSpatialFrequencyRange(2) <= 70)
        xTicks = [0.01 0.03 0.1 0.3 1 3 10 30 60];
    else
        xTicks = [0.01 0.03 0.1 0.3 1 3 10 30 60 100];
    end
    
    visualizedSpatialFrequencyRange = [visualizedSpatialFrequencyRange(1) visualizedSpatialFrequencyRange(2)];
    
    grid(ax, 'on'); box(ax, 'off');
    set(ax, 'XLim', visualizedSpatialFrequencyRange, 'YLim', [1*ff.axisOffsetFactor 0.5]);
    
    set(ax, 'XTick', xTicks, 'YTick', 0:0.1:1, ...
            'XTickLabel', xTicks, 'YTickLabel', 0:0.1:1);
    set(ax, 'TickDir', 'both')
    set(ax, 'XScale', 'log');
    

    % xy axis labels
    xtickangle(ax, 0);
    ytickangle(ax, 0);
    if (~noXLabel)
        xlabel(ax, 'spatial frequency (c/deg)', 'FontAngle', ff.axisFontAngle);
    else
        set(ax, 'XTickLabel', {});
    end
    if (~noYLabel)
        ylabel(ax, 'STF magnitude', 'FontAngle', ff.axisFontAngle);
    end
    
    
    if (~noYTickLabel)
        set(ax, 'YTick', 0:0.1:1, 'YTickLabel', sprintf('%1.1f\n', 0:0.1:1));
    else
        set(ax, 'YTick', []);
    end

    if (~noXTickLabel)
    else
        set(ax, 'XTickLabel', {});
    end


    % Font size
    set(ax, 'FontSize', ff.fontSize);

    % axis color and width
    set(ax, 'XColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    set(ax, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);

    box(ax, 'off');

    yyaxis(ax, 'left')
    %hold(ax, 'on')

    % The difference between the 2 cone-isolating STF phase spectra
    theConeIsolatingPhaseSpectrumDifferential = unwrapPhaseSpectum(squeeze(theSTFphases(:,2)), squeeze(theSTFphases(:,3)));

    plot(ax, sfSupportCPD, theConeIsolatingPhaseSpectrumDifferential, '-', 'LineWidth', ff.lineWidth, ...
            'MarkerEdgeColor', [0 0 0.8], ...
            'MarkerFaceColor', [0.5 0.5 1], ...
            'Color', [0 0.3 0.6], ...
            'MarkerSize', ff.markerSize/2, ...
            'LineWidth', ff.lineWidth*2);
    hold(ax, 'on');
    plot(ax, sfSupportCPD, theConeIsolatingPhaseSpectrumDifferential, '-', 'LineWidth', ff.lineWidth, ...
            'MarkerEdgeColor', [0 0 0.8], ...
            'MarkerFaceColor', [0.5 0.5 1], ...
            'Color', [0 0.7 0.9], ...
            'MarkerSize', ff.markerSize/2, ...
            'LineWidth', ff.lineWidth);
    hold(ax, 'off')
    minPhase = min(theConeIsolatingPhaseSpectrumDifferential);
    maxPhase = max(theConeIsolatingPhaseSpectrumDifferential);

    phaseRange(1) = (sign(minPhase)*round(abs(minPhase)/30)-1) * 30;
    phaseRange(2) = phaseRange(1) + 180 + 60;

    while (phaseRange(2) < maxPhase)
       phaseRange = phaseRange + 30; 
    end

    phaseTicks = (phaseRange(1)+30):30:(phaseRange(2)-30);
    phaseTickLabels = sprintf('%2.0f\n', phaseTicks);

    set(ax, 'YLim', phaseRange, 'YTick', phaseTicks, 'YTickLabel', phaseTickLabels, 'YColor', [0 0 0.8]);
    ylabel(ax, 'STF differential phase (L-M)', 'FontAngle', ff.axisFontAngle);

    set(ax, 'YColor', [0 0.6 0.8], 'LineWidth', ff.axisLineWidth);

    % plot title
    if (~isempty(plotTitle))
        title(ax, plotTitle, 'Color', ff.titleColor, 'FontSize', ff.titleFontSize, 'FontWeight', ff.titleFontWeight);
    end

    % legends
    yyaxis(ax, 'right')
    if (~isempty(theLegends))
        legend(ax, pH, theLegends, ...
            'Location', 'EastOutside', 'NumColumns', 1, ...
            'FontSize', ff.legendFontSize, 'Box', 'off');
    end

end

function theUnwrappedPhaseSpectrum = unwrapPhaseSpectum(theLconePhaseSpectrum, theMconePhaseSpectrum)
    theLMdifferentialPhaseSpectrum = theLconePhaseSpectrum - theMconePhaseSpectrum;
    theUnwrappedPhaseSpectrum = unwrap(theLMdifferentialPhaseSpectrum/180*pi)/pi*180;
end

function p = shadedAreaBetweenTwoLines(ax,x, y1, y2, ...
     faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)

    y1 = reshape(y1, size(x));
    y2 = reshape(y2, size(x));

    x = [x  x(end)  fliplr(x)  x(1)];
    y = [y1 y2(end) fliplr(y2) y2(1)];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    p = patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, ...
        'FaceAlpha', faceAlpha, 'LineWidth', lineWidth, 'LineStyle', lineStyle);
 end