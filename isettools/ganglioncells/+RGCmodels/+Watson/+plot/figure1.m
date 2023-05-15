% Generate Figure 1 of the Watson (2014) paper, which plots the variation
% in cone density with eccentricity along the 4 principal meridians.
function figure1()
    % Define ecc range in degrees
    eccDegs = logspace(log10(0.01), log10(90), 50);
    
    % Retrieve all meridians
    examinedMeridians = RGCmodels.Watson.constants.indexedMeridians;
    
    % Compute cone density along each meridian
    useParfor = true;
    coneDensityDegs2 = zeros(numel(examinedMeridians), numel(eccDegs));
    for k = 1:numel(examinedMeridians)
        [~,~, coneDensityDegs2(k,:)] = ...
            RGCmodels.Watson.compute.coneSpacingAlongMeridianInRightEyeVisualField(eccDegs, examinedMeridians{k}, useParfor);
    end
    
    % Generate figure
    hFig = figure(1); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);

    theMeridianLabels = cell(1,numel(examinedMeridians));
    for k = 1:numel(examinedMeridians)
        % Retrieve color for meridian
        meridianColor = RGCmodels.Watson.plot.colorForMeridian(examinedMeridians{k}); 
        theMeridianLabels{k} = examinedMeridians{k};
        plot(theAxes{1,1}, eccDegs, squeeze(coneDensityDegs2(k,:)), 'k-', 'Color', meridianColor, 'LineWidth', 1.5);
        hold(theAxes{1,1}, 'on');
    end

    % Plot the max cone density as a disk
    plot(theAxes{1,1},0.01, RGCmodels.Watson.constants.dc0, 'ko', 'MarkerSize', 14, 'LineWidth', 1.0, 'MarkerFaceColor', 0.8*[1 1 1]);
    
    % Label the meridians
    legend(theAxes{1,1},theMeridianLabels,  'Location', 'SouthWest');

    % Axes and limits
    set(theAxes{1,1}, 'XLim', [0.01*(1+3*ff.axisOffsetFactor) 100], 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100]);
    set(theAxes{1,1}, 'YLim', [100*(1+3*ff.axisOffsetFactor) 20000], 'YScale', 'log', 'YTick', [100 300 1000 3000 10000]);
    
    xlabel(theAxes{1,1}, 'eccentricity (degs)');
    ylabel(theAxes{1,1},'density (cones/deg^2)');

    % Grid
    grid(theAxes{1,1}, 'on'); box(theAxes{1,1}, 'off');

    % Ticks
    set(theAxes{1,1}, 'TickDir', 'both');

    % Font size
    set(theAxes{1,1}, 'FontSize', ff.fontSize+4);

    % axis color and width
    set(theAxes{1,1}, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
end



