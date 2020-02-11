function hFig = generateFigure14(obj, hFig)
% Generate Figure 14 of the Watson 2014 paper
%
% Syntax:
%   WatsonRGCCalc = WatsonRGCModel();
%   WatsonRGCCalc.generateFigure14();
%
% Description:
%   Generate Figure 14 of Watson (2014) which plots the ratio of midget to
%   cones as a function of eccentricity for all 4 meridians.
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

    figureNumber = '14';
    figure(hFig); clf;
    set(hFig, 'Color', obj.figurePrefs.backgroundColor, 'Name', sprintf('Figure %s of %s', figureNumber, obj.paperTitleShort));
    
    meridianName = 'temporal meridian';
    midgetRGCRFDensity = obj.midgetRGCRFDensity(obj.eccDegs, meridianName, 'RFs per deg2');
    [~, coneRFDensity] = obj.coneRFSpacingAndDensity(obj.eccDegs, meridianName, 'Cones per deg2');
    plot(obj.eccDegs, midgetRGCRFDensity./coneRFDensity, 'r-', 'LineWidth', obj.figurePrefs.lineWidth); hold on
        
    meridianName = 'superior meridian';
    midgetRGCRFDensity = obj.midgetRGCRFDensity(obj.eccDegs, meridianName, 'RFs per deg2');
    [~, coneRFDensity] = obj.coneRFSpacingAndDensity(obj.eccDegs, meridianName, 'Cones per deg2');
    plot(obj.eccDegs, midgetRGCRFDensity./coneRFDensity, 'b-', 'LineWidth', obj.figurePrefs.lineWidth);
    
    meridianName = 'nasal meridian';
    midgetRGCRFDensity = obj.midgetRGCRFDensity(obj.eccDegs, meridianName, 'RFs per deg2');
    [~, coneRFDensity] = obj.coneRFSpacingAndDensity(obj.eccDegs, meridianName, 'Cones per deg2');
    plot(obj.eccDegs, midgetRGCRFDensity./coneRFDensity, 'g-', 'LineWidth', obj.figurePrefs.lineWidth);
    
    meridianName = 'inferior meridian';
    midgetRGCRFDensity = obj.midgetRGCRFDensity(obj.eccDegs, meridianName, 'RFs per deg2');
    [~, coneRFDensity] = obj.coneRFSpacingAndDensity(obj.eccDegs, meridianName, 'Cones per deg2');
    plot(obj.eccDegs, midgetRGCRFDensity./coneRFDensity, 'k-', 'LineWidth', obj.figurePrefs.lineWidth);
    
    legend({'temporal', 'superior', 'nasal', 'inferior'}, 'Location', 'SouthWest');

    xlabel('eccentricity (degs)', 'FontAngle', obj.figurePrefs.fontAngle);
    ylabel('midget RGC / cone ratio', 'FontAngle', obj.figurePrefs.fontAngle);
    set(gca, 'XLim', [0.01 80], 'YLim', [0.02 2.2], ...
        'XScale', 'log', 'YScale', 'log', ...
        'XTick', [0.1 0.5 1 5 10 50 100], 'YTick', [0.05 0.1 0.2 0.5 1 2], ...
        'FontSize', obj.figurePrefs.fontSize);
    grid(gca, obj.figurePrefs.grid);
    
    title('Ratio of midget RGC RFs to cones as a function of eccentricity');
end