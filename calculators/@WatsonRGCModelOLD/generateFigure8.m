function hFig = generateFigure8(obj, hFig)
% Generate Figure 8 of the Watson 2014 paper
%
% Syntax:
%   WatsonRGCCalc = WatsonRGCModel();
%   WatsonRGCCalc.generateFigure8();
%
% Description:
%   Generate Figure 8 of Watson (2014) which plots the fraction of midget to
%   total RGC RF density as a function of eccentricity.
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

    figureNumber = '8';
    figure(hFig); clf;
    set(hFig, 'Color', obj.figurePrefs.backgroundColor, 'Name', sprintf('Figure %s of %s', figureNumber, obj.paperTitleShort));
    midgetRGCFraction = obj.midgetRGCFraction(obj.eccDegs);
    plot(obj.eccDegs, midgetRGCFraction, 'r-', 'LineWidth', obj.figurePrefs.lineWidth);
    xlabel('eccentricity (degs)', 'FontAngle', obj.figurePrefs.fontAngle);
    ylabel('midget RGC RF fraction', 'FontAngle', obj.figurePrefs.fontAngle);
    set(gca, 'XLim', [0 100], 'YLim', [0 1], ...
        'XTick', 0:5:100, 'YTick', [0:0.1:1.0], ...
        'FontSize', obj.figurePrefs.fontSize);
    grid(gca, obj.figurePrefs.grid);
    
    title('Fraction of midget RGC RFs as a function of eccentricity');
end