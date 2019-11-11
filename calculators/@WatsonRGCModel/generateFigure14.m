function hFig = generateFigure14(obj)
% Generate Figure 14 of the Watson 2014 paper
%
% Syntax:
%   WatsonRGCCalc = WatsonRGCModel();
%   WatsonRGCCalc.generateFigure14();
%
% Description:
%   Generate Figure 14 of Watson (2014) which plots the ratio of midget to
%   cones as a function of eccentricity.
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
    hFig = figure(); clf;
    set(hFig, 'Color', obj.figurePrefs.backgroundColor, 'Name', sprintf('Figure %s of %s', figureNumber, obj.paperTitleShort));
    eccDegs = 0.1:0.002:70;
    
    meridianName = 'temporal meridian';
    midgetRGCRFDensity = obj.midgetRGCRFDensity(eccDegs, meridianName, 'RFs per mm2');
    [~, coneRFDensity] = obj.coneRFSpacingAndDensity(eccDegs, meridianName, 'mm');
    plot(eccDegs, midgetRGCRFDensity./coneRFDensity, 'r-', 'LineWidth', obj.figurePrefs.lineWidth); hold on
        
    meridianName = 'superior meridian';
    midgetRGCRFDensity = obj.midgetRGCRFDensity(eccDegs, meridianName, 'RFs per mm2');
    [~, coneRFDensity] = obj.coneRFSpacingAndDensity(eccDegs, meridianName, 'mm');
    plot(eccDegs, midgetRGCRFDensity./coneRFDensity, 'b-', 'LineWidth', obj.figurePrefs.lineWidth);
    
    meridianName = 'nasal meridian';
    midgetRGCRFDensity = obj.midgetRGCRFDensity(eccDegs, meridianName, 'RFs per mm2');
    [~, coneRFDensity] = obj.coneRFSpacingAndDensity(eccDegs, meridianName, 'mm');
    plot(eccDegs, midgetRGCRFDensity./coneRFDensity, 'g-', 'LineWidth', obj.figurePrefs.lineWidth);
    
    meridianName = 'inferior meridian';
    midgetRGCRFDensity = obj.midgetRGCRFDensity(eccDegs, meridianName, 'RFs per mm2');
    [~, coneRFDensity] = obj.coneRFSpacingAndDensity(eccDegs, meridianName, 'mm');
    plot(eccDegs, midgetRGCRFDensity./coneRFDensity, 'k-', 'LineWidth', obj.figurePrefs.lineWidth);
    
    legend({'temporal', 'superior', 'nasal', 'inferior'}, 'Location', 'SouthWest');

    xlabel('eccentricity (degs)', 'FontAngle', obj.figurePrefs.fontAngle);
    ylabel('midget RGC / cone ratio', 'FontAngle', obj.figurePrefs.fontAngle);
    set(gca, 'XLim', [0.1 80], 'YLim', [0.02 2.2], ...
        'XScale', 'log', 'YScale', 'log', ...
        'XTick', [0.5 1 5 10 50], 'YTick', [0.05 0.1 0.2 0.5 1 2], ...
        'FontSize', obj.figurePrefs.fontSize);
    grid(gca, obj.figurePrefs.grid);
    
    title('Fraction of midget RGC RFs as a function of eccentricity');
end