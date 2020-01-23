function hFig = generateFigure1(obj, hFig)
% Generate Figure 1 of the Watson 2014 paper
%
% Syntax:
%   WatsonRGCCalc = WatsonRGCModel();
%   WatsonRGCCalc.generateFigure1();
%
% Description:
%   Generate Figure 14 of Watson (2014) which plots cone density as a 
%   function of eccentricity for all 4 meridians.
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

    figureNumber = '1';
    figure(hFig); clf;
    set(hFig, 'Color', obj.figurePrefs.backgroundColor, 'Name', sprintf('Figure %s of %s', figureNumber, obj.paperTitleShort));
    
    
    meridianName = 'temporal meridian';
    [~, coneRFDensity] = obj.coneRFSpacingAndDensity(obj.eccDegs, meridianName, 'Cones per deg2');
    plot(obj.eccDegs, coneRFDensity, 'r-', 'LineWidth', obj.figurePrefs.lineWidth); hold on
        
    meridianName = 'superior meridian';
    [~, coneRFDensity] = obj.coneRFSpacingAndDensity(obj.eccDegs, meridianName, 'Cones per deg2');
    plot(obj.eccDegs, coneRFDensity, 'b-', 'LineWidth', obj.figurePrefs.lineWidth);
    
    meridianName = 'nasal meridian';
    [~, coneRFDensity] = obj.coneRFSpacingAndDensity(obj.eccDegs, meridianName, 'Cones per deg2');
    plot(obj.eccDegs, coneRFDensity, 'g-', 'LineWidth', obj.figurePrefs.lineWidth); 
    
    meridianName = 'inferior meridian';
    [~, coneRFDensity] = obj.coneRFSpacingAndDensity(obj.eccDegs, meridianName, 'Cones per deg2');
    plot(obj.eccDegs, coneRFDensity, 'k-', 'LineWidth', obj.figurePrefs.lineWidth);
    
    plot([obj.eccDegs(2) obj.eccDegs(end)], obj.dc0*[1 1], 'k--', 'LineWidth', 1.5);
    legend({'temporal', 'superior', 'nasal', 'inferior', '(0,0)'}, 'Location', 'SouthWest');

    xlabel('eccentricity (degs)', 'FontAngle', obj.figurePrefs.fontAngle);
    ylabel('density (cones/deg^2)', 'FontAngle', obj.figurePrefs.fontAngle);
    set(gca, 'XLim', [0.005 100], 'YLim', [100 20000], ...
        'XScale', 'log', 'YScale', 'log', ...
        'XTick', [0.1 0.5 1 5 10 50 100], ...
        'YTick', [100 1000 10000], 'YTickLabel', {'100', '1000', '10000'},...
        'FontSize', obj.figurePrefs.fontSize);
    grid(gca, obj.figurePrefs.grid);
    
    title('Cone density as a function of eccentricity (Curcio et al, 1990)');
    drawnow;
end
