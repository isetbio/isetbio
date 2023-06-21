% Generate Figure 14 of the Watson (2014) paper, which plots the ratio of
% midget RGCs: cones as a function of eccentricity along the 4 principal meridians.
function hFig = figure14()
    % Define ecc range in degrees
    eccDegs = logspace(log10(0.01), log10(90), 50);
    
    % Retrieve all meridians
    examinedMeridians = RGCmodels.Watson.constants.indexedMeridians;
    
    % Compute cone density along each meridian
    midgetRGCtoConeRatio = zeros(numel(examinedMeridians), numel(eccDegs));
    useParfor = true;
    for k = 1:numel(examinedMeridians)
         [~,~, coneDensityDegs2] = ...
            RGCmodels.Watson.compute.coneSpacingAlongMeridianInRightEyeVisualField(eccDegs, examinedMeridians{k}, useParfor);
        [~,~, ONandOFFmidgetRGCRFDensityDegs2] = ...
            RGCmodels.Watson.compute.midgetRGCRFSpacingAlongMeridianInRightEyeVisualField(eccDegs, examinedMeridians{k});
        midgetRGCtoConeRatio(k,:) = ONandOFFmidgetRGCRFDensityDegs2 ./ coneDensityDegs2;
    end
    
    % Only one of subtype (half of ON+OFF)
    midgetRGCRFDensityDegs2 = ONandOFFmidgetRGCRFDensityDegs2 / 2;
    midgetRGCtoConeRatio =  midgetRGCtoConeRatio /2;

    % Generate figure
    hFig = figure(14); clf;
    ff = MSreadyPlot.figureFormat('1x1 small wide');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);

    theMeridianLabels = cell(1,numel(examinedMeridians));
    for k = 1:numel(examinedMeridians)
        % Retrieve color for meridian
        meridianColor = RGCmodels.Watson.plot.colorForMeridian(examinedMeridians{k});
        theMeridianLabels{k} = examinedMeridians{k};
        plot(theAxes{1,1},eccDegs, squeeze(midgetRGCtoConeRatio(k,:)), 'k-', 'Color', meridianColor, 'LineWidth', ff.lineWidth);
        hold(theAxes{1,1}, 'on');
    end

    % Plot the max cone density as a disk
    plot(theAxes{1,1},0.01, max(midgetRGCtoConeRatio(:)), 'ko', ...
        'MarkerSize', ff.markerSize, 'LineWidth', 1.0, 'MarkerFaceColor', 0.8*[1 1 1]);
  

    % Label the meridians
    legend(theAxes{1,1},theMeridianLabels, 'Location', 'SouthWest', 'FontSize', ff.legendFontSize);
    legend boxoff

    % Axes and limits
    set(theAxes{1,1}, 'XLim', [0.01*(1+3*ff.axisOffsetFactor) 100], 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100]);
    set(theAxes{1,1}, 'YLim', [0.01*(1+3*ff.axisOffsetFactor) 1.1], 'YScale', 'log', ...
        'YTick', [0.01 0.02 0.05 0.1 0.2 0.5 1 2], ...
        'YTickLabel', {'.01', '.02', '.05', '.1', '.2', '.5', '1', '2'});
    xtickangle(theAxes{1,1}, 0);

    xlabel(theAxes{1,1},'eccentricity (degs)', 'FontAngle', ff.axisFontAngle);
    ylabel(theAxes{1,1}, 'midget RGC RF / cone ratio', 'FontAngle', ff.axisFontAngle);

    % Grid
    grid(theAxes{1,1}, 'on'); box(theAxes{1,1}, 'off');

    % Ticks
    set(theAxes{1,1}, 'TickDir', 'both');

    % Font size
    set(theAxes{1,1}, 'FontSize', ff.fontSize);

    % axis color and width
    set(theAxes{1,1}, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
end



