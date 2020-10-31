% Generate Figure 14 of the Watson (2014) paper, which plots the ratio of
% midget RGCs: cones as a function of eccentricity along the 4 principal meridians.
function figure14()
    % Define ecc range in degrees
    eccDegs = logspace(log10(0.01), log10(80), 50);
    
    % Retrieve all meridians
    examinedMeridians = RGCmodels.Watson.constants.indexedMeridians;
    
    % Compute cone density along each meridian
    midgetRGCtoConeRatio = zeros(numel(examinedMeridians), numel(eccDegs));
    for k = 1:numel(examinedMeridians)
         [~,~, coneDensityDegs2] = ...
            RGCmodels.Watson.compute.coneSpacingAlongMeridianInRightEyeVisualField(eccDegs, examinedMeridians{k});
        [~,~, midgetRGCRFDensityDegs2] = ...
            RGCmodels.Watson.compute.midgetRGCRFSpacingAlongMeridianInRightEyeVisualField(eccDegs, examinedMeridians{k});
        midgetRGCtoConeRatio(k,:) = midgetRGCRFDensityDegs2 ./ coneDensityDegs2;
    end
    
    % Generate figure
    hFig = figure(14); clf;
    set(hFig, 'Position', [10 10 1000 500], 'Color', [1 1 1]);
    theMeridianLabels = cell(1,numel(examinedMeridians));
    for k = 1:numel(examinedMeridians)
        % Retrieve color for meridian
        meridianColor = RGCmodels.Watson.plot.colorForMeridian(examinedMeridians{k}); hold on
        theMeridianLabels{k} = examinedMeridians{k};
        plot(eccDegs, squeeze(midgetRGCtoConeRatio(k,:)), 'k-', 'Color', meridianColor, 'LineWidth', 1.5);
    end

    % Label the meridians
    legend(theMeridianLabels);
    % Axes and limits
    set(gca, 'XLim', [0.01 100], 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100]);
    set(gca, 'YLim', [0.02 3], 'YScale', 'log', 'YTick', [0.05 0.1 0.2 0.5 1 2]);
    set(gca, 'FontSize', 20);
    xlabel('eccentricity (degs)');
    ylabel('midget RGC RF / cone ratio');
    grid on;
end



