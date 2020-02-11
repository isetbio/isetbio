function hFig = generateFigureA1(obj, hFig)
% Generate Figure A1 of the Watson 2014 paper
%
% Syntax:
%   WatsonRGCCalc = WatsonRGCModel();
%   WatsonRGCCalc.generateFigureA1();
%
% Description:
%   Generate Figure A1 of Watson (2014) which plots the relationhip between retinal distance from the optic axis in mm and degs as a
%   function of eccentricity
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

    figureNumber = 'A1-left';
    figure(hFig); clf;
    set(hFig, 'Color', obj.figurePrefs.backgroundColor, 'Name', sprintf('Figure %s of %s', figureNumber, obj.paperTitleShort));
    eccMM = obj.rhoDegsToMMs(obj.eccDegs);
    linearApproximationDegsToMMs = 268/1000;
    plot(obj.eccDegs, eccMM, 'k-', 'LineWidth', obj.figurePrefs.lineWidth); hold on
    plot(obj.eccDegs, linearApproximationDegsToMMs*obj.eccDegs, 'r-', 'LineWidth', obj.figurePrefs.lineWidth);
    legend({'deg->mm (ecc)', 'linear approximation'}, 'Location', 'SouthEast');
    axis 'square';
    xlabel('eccentricity, r'' (degs)', 'FontAngle', obj.figurePrefs.fontAngle);
    ylabel('eccentricity, r'' (mm)', 'FontAngle', obj.figurePrefs.fontAngle);
    set(gca, 'XLim', [0 100], 'YLim', [0 27], 'XTick', 0:20:100, 'YTick', 0:5:30, 'FontSize', obj.figurePrefs.fontSize);
    grid(gca, obj.figurePrefs.grid);
    title('Relationhip between retinal distance from the optic axis in mm and degs');
    
%     figureNumber = 'A1-right';
%     hFig = figure(); clf;
%     set(hFig, 'Color', obj.figurePrefs.backgroundColor, 'Name', sprintf('Figure %s of %s', figureNumber, obj.paperTitleShort));
%     eccMM = 0:0.1:22;
%     eccDegs = obj.rhoMMsToDegs(eccMM);
%     plot(eccMM, eccDegs, 'k-', 'LineWidth', obj.figurePrefs.lineWidth); hold on
%     plot(eccMM, 1 / linearApproximationDegsToMMs * eccMM, 'r-', 'LineWidth', obj.figurePrefs.lineWidth);
%     legend({'mm->deg (ecc)', 'linear approximation'}, 'Location', 'SouthEast');
%     axis 'square';
%     xlabel('eccentricity, r'' (mm)', 'FontAngle', obj.figurePrefs.fontAngle);
%     ylabel('eccentricity, r'' (degs)', 'FontAngle', obj.figurePrefs.fontAngle);
%     set(gca, 'YLim', [0 100], 'XLim', [0 22], 'YTick', 0:20:100, 'XTick', 0:5:30, 'FontSize', obj.figurePrefs.fontSize);
%     grid(gca, obj.figurePrefs.grid);
%     title('Relationhip between retinal distance from the optic axis in mm and degs');
    
end