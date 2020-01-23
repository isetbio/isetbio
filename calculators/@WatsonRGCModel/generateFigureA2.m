function hFig = generateFigureA2(obj, hFig)
% Generate Figure A2 of the Watson 2014 paper
%
% Syntax:
%   WatsonRGCCalc = WatsonRGCModel();
%   WatsonRGCCalc.generateFigureA2();
%
% Description:
%   Generate Figure A2 of Watson (2014) which plots the ratio of area in
%   mm2 to area in deg^2 as a function of eccentricity.
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

    figureNumber = 'A2';
    figure(hFig); clf;
    set(hFig, 'Color', obj.figurePrefs.backgroundColor, 'Name', sprintf('Figure %s of %s', figureNumber, obj.paperTitleShort));
    plot(obj.eccDegs, obj.alpha(eccDegs), 'r-', 'LineWidth', obj.figurePrefs.lineWidth);
    xlabel('eccentricity, r'' (degs)', 'FontAngle', obj.figurePrefs.fontAngle);
    ylabel('area ratio, \alpha (mm^2/deg^2)', 'FontAngle', obj.figurePrefs.fontAngle);
    set(gca, 'XLim', [0 100], 'YLim', [0.01 0.08], 'FontSize', obj.figurePrefs.fontSize);
    grid(gca, obj.figurePrefs.grid);
    title('Area ratio as a function of angular distance from the optic axis');
end