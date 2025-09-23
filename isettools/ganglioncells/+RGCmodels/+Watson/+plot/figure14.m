% Generate Figure 14 of the Watson (2014) paper, which plots the ratio of
% midget RGCs: cones as a function of eccentricity along the 4 principal meridians.
function hFig = figure14(varargin)

    % Parse input
    p = inputParser;
    p.addParameter('inverseRatio', false, @islogical);
    p.parse(varargin{:});

    inverseRatio = p.Results.inverseRatio;

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
    set(hFig, 'Color', [1 1 1]);
    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);

    theRatio = midgetRGCtoConeRatio;
    if (inverseRatio)
         theRatio = 1./ theRatio;
    end

    theMeridianLabels = cell(1,numel(examinedMeridians));
    pHandles = zeros(1,numel(examinedMeridians));
    for k = 1:numel(examinedMeridians)
        % Retrieve color for meridian
        meridianColor = RGCmodels.Watson.plot.colorForMeridian(examinedMeridians{k});
        theMeridianLabels{k} = examinedMeridians{k};
        pHandles(k) = plot(theAxes{1,1},eccDegs, squeeze(theRatio(k,:)), 'k-', 'Color', meridianColor, 'LineWidth', ff.lineWidth);
        hold(theAxes{1,1}, 'on');
    end

    % Axes scaling and ticks
    set(theAxes{1,1}, 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100]);
    
    % Finalize figure using the Publication-Ready format
    XLims = [0.01 100];
    
    if (inverseRatio)

        % Plot the min cone density as a disk
        plot(theAxes{1,1},0.01, min(theRatio(:)), 'ko', ...
            'MarkerSize', ff.markerSize, 'LineWidth', 1.0, 'MarkerFaceColor', 0.8*[1 1 1]);

        % Label the meridians
        legend(theAxes{1,1}, pHandles, theMeridianLabels, 'Location', 'NorthWest');

        set(theAxes{1,1}, 'YScale', 'log', ...
            'YTick', [0.5 1 1.7 3 5 10 18 30 55 100 200], ...
            'YTickLabel', {'.5', '1', '1.7', '3', '5', '10', '18', '30', '55', '100', '200'});
    
        YLims = [1 200];
        PublicationReadyPlotLib.offsetAxes(theAxes{1,1},ff, XLims, YLims);
        PublicationReadyPlotLib.labelAxes(theAxes{1,1},ff, 'eccentricity (degs)', 'cone / midget RGC RF ratio');

    else

        % Plot the max cone density as a disk
        plot(theAxes{1,1},0.01, max(theRatio(:)), 'ko', ...
            'MarkerSize', ff.markerSize, 'LineWidth', 1.0, 'MarkerFaceColor', 0.8*[1 1 1]);

        % Label the meridians
        legend(theAxes{1,1}, pHandles, theMeridianLabels, 'Location', 'NorthWest');

        set(theAxes{1,1}, 'YScale', 'log', ...
            'YTick', [0.01 0.02 0.05 0.1 0.2 0.5 1 2], ...
            'YTickLabel', {'.01', '.02', '.05', '.1', '.2', '.5', '1', '2'});
    
        YLims = [0.01 1.1];
        PublicationReadyPlotLib.offsetAxes(theAxes{1,1},ff, XLims, YLims);
        PublicationReadyPlotLib.labelAxes(theAxes{1,1},ff, 'eccentricity (degs)', 'midget RGC RF / cone ratio');

    end

    PublicationReadyPlotLib.applyFormat(theAxes{1,1},ff);
end



