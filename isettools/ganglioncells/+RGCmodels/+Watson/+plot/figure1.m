% Generate Figure 1 of the Watson (2014) paper, which plots the variation
% in cone density with eccentricity along the 4 principal meridians.
function figure1()
    % Define ecc range in degrees
    eccDegs = logspace(log10(0.01), log10(90), 50);
    
    % Retrieve all meridians
    examinedMeridians = RGCmodels.Watson.constants.indexedMeridians;
    
    % Compute cone density along each meridian
    coneDensityDegs2 = zeros(numel(examinedMeridians), numel(eccDegs));
    for k = 1:numel(examinedMeridians)
        [~,~, coneDensityDegs2(k,:)] = ...
            RGCmodels.Watson.compute.coneSpacingAlongMeridianInRightEyeVisualField(eccDegs, examinedMeridians{k});
    end
    
    % Generate figure
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1000 500], 'Color', [1 1 1]);
    theMeridianLabels = cell(1,numel(examinedMeridians));
    for k = 1:numel(examinedMeridians)
        % Retrieve color for meridian
        meridianColor = RGCmodels.Watson.plot.colorForMeridian(examinedMeridians{k}); hold on
        theMeridianLabels{k} = examinedMeridians{k};
        plot(eccDegs, squeeze(coneDensityDegs2(k,:)), 'k-', 'Color', meridianColor, 'LineWidth', 1.5);
    end
    % Plot the max cone density as a disk
    plot(0.01, RGCmodels.Watson.constants.dc0, 'ko', 'MarkerSize', 14, 'MarkerFaceColor', [1 1 1]);
    % Label the meridians
    legend(theMeridianLabels);
    % Axes and limits
    set(gca, 'XLim', [0.01 100], 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100]);
    set(gca, 'YLim', [100 20000], 'YScale', 'log', 'YTick', [100 300 1000 3000 10000]);
    set(gca, 'FontSize', 20);
    xlabel('eccentricity (degs)');
    ylabel('density (cones/deg^2)');
    grid on;
end



