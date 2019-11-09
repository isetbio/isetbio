function hFig = generateFigure5(obj)
% Generate Figure 5 of the Watson 2014 paper
%
% Syntax:
%   WatsonRGCCalc = WatsonRGCModel();
%   WatsonRGCCalc.generateFigure5();
%
% Description:
%   Generate Figure 5 of Watson (2014) which plots the total RGC RF density
%   as a function of eccentricity for all 4 meridians.
%
% Inputs:
%    obj                       - The WatsonRGCModel object
% Outputs:
%    hFig                      - The generated figure handle
%
% References:
%    Watson (2014). 'A formula for human RGC receptive field density as
%    a function of visual field location', JOV (2014), 14(7), 1-17.
%
% History:
%    11/8/19  NPC, ISETBIO Team     Wrote it.

    figureNumber = '5';
    hFig = figure(); clf;
    set(hFig, 'Color', obj.figurePrefs.backgroundColor, 'Name', sprintf('Figure %s of %s', figureNumber, obj.paperTitleShort));
    eccDegs = 0.1:0.1:100;
    % Plot the total RGC RF density along the temporal meridian
    plot(eccDegs, obj.totalRGCRFDensity(eccDegs, 'temporal meridian', 'RFs per deg2'), ...
        'r-', 'LineWidth', obj.figurePrefs.lineWidth); hold on;
    plot(eccDegs, obj.totalRGCRFDensity(eccDegs, 'superior meridian', 'RFs per deg2'), ...
        'b-', 'LineWidth', obj.figurePrefs.lineWidth);
    plot(eccDegs, obj.totalRGCRFDensity(eccDegs, 'nasal meridian', 'RFs per deg2'), ...
        'g-', 'LineWidth', obj.figurePrefs.lineWidth);
    plot(eccDegs, obj.totalRGCRFDensity(eccDegs, 'inferior meridian', 'RFs per deg2'), ...
        'k-', 'LineWidth', obj.figurePrefs.lineWidth);
    legend({'temporal', 'superior', 'nasal', 'inferior'}, 'Location', 'SouthWest');
    xlabel('eccentricity (degs)', 'FontAngle', obj.figurePrefs.fontAngle);
    ylabel('total RGC density (# of RFs/deg^2)', 'FontAngle', obj.figurePrefs.fontAngle);
    set(gca, 'XLim', [0.07 120], 'YLim', [3 35*1000], ...
        'XTick', [0.1 0.5 1 5 10 50 100], 'YTick', [10 100 1000 10000], ...
        'XScale', 'log', 'YScale', 'log', 'FontSize', obj.figurePrefs.fontSize);
    grid(gca, obj.figurePrefs.grid);
    title('Total RGC RF density as a function of eccentricity for four meridians');
end