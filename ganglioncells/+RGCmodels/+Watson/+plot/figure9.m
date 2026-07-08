% Generate Figure 9 of the Watson (2014) paper, which plots the variation
% in midget RGC RF with eccentricity along the 4 principal meridians.
function hFig = figure9()
    % Define ecc range in degrees
    eccDegs = logspace(log10(0.01), log10(90), 50);
    
    % Retrieve all meridians
    examinedMeridians = RGCmodels.Watson.constants.indexedMeridians;
    
    % Compute cone density along each meridian
    midgetRGCRFDensityDegs2 = zeros(numel(examinedMeridians), numel(eccDegs));
    for k = 1:numel(examinedMeridians)
        [~,~, ONandOFFmidgetRGCRFDensityDegs2(k,:)] = ...
            RGCmodels.Watson.compute.midgetRGCRFSpacingAlongMeridianInRightEyeVisualField(eccDegs, examinedMeridians{k});
    end
    
    % Only one of subtype (half of ON+OFF)
    midgetRGCRFDensityDegs2 = ONandOFFmidgetRGCRFDensityDegs2 / 2;

    % Generate figure
    hFig = figure(9); clf;
    set(hFig, 'Color', [1 1 1]);
    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);

    theMeridianLabels = cell(1,numel(examinedMeridians));
    for k = 1:numel(examinedMeridians)
        % Retrieve color for meridian
        meridianColor = RGCmodels.Watson.plot.colorForMeridian(examinedMeridians{k}); hold on
        theMeridianLabels{k} = examinedMeridians{k};
        plot(theAxes{1,1}, eccDegs, squeeze(midgetRGCRFDensityDegs2(k,:)), 'k-', 'Color', meridianColor, 'LineWidth', ff.lineWidth);
    end

    % Plot the max midget RGC density as a disk
    plot(theAxes{1,1}, 0.01, RGCmodels.Watson.constants.dc0, 'ko', ...
        'MarkerSize', ff.markerSize, 'LineWidth', 1.0, 'MarkerFaceColor', 0.8*[1 1 1]);

    % Label the meridians
    legend(theAxes{1,1}, theMeridianLabels,  'Location', 'SouthWest');

    % Axes scaling and ticks
    set(theAxes{1,1}, 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100]);
    set(theAxes{1,1}, 'YScale', 'log', ...
        'YTick', [1 10 100 1000 10000], ...
        'YTickLabel', {'1', '10', '100', '1k', '10k'});

    % Finalize figure using the Publication-Ready format
    XLims = [0.01 100]; YLims = [1 20000];
    PublicationReadyPlotLib.offsetAxes(theAxes{1,1},ff, XLims, YLims);
    PublicationReadyPlotLib.labelAxes(theAxes{1,1},ff, 'eccentricity (degs)', 'density (midget RGC RFs/deg^2)');
    PublicationReadyPlotLib.applyFormat(theAxes{1,1},ff);
end



