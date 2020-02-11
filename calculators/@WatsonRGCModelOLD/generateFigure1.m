function hFig = generateFigure1(obj, hFig, varargin)
% Generate Figure 1 of the Watson 2014 paper
%
% Syntax:
%   WatsonRGCCalc = WatsonRGCModel();
%   WatsonRGCCalc.generateFigure1();
%
% Description:
%   Generate Figure 1 of Watson (2014) which plots cone density in cones/deg^2 
%   as a function of eccentricity for all 4 meridians.
%
% Inputs:
%    obj                            - The WatsonRGCModel object
%    hFig                           - Figure handle
%    
% Optional Inputs:
%    whichEye                       - Which eye, 'left' or 'right'. The default 
%                                     is 'right' to match Watson (2014).
%
%    eccentricityInMMInsteadOfDegs  - Logical. Whether to plot as a function
%                                     of retinal distance specified in mm. 
%                                     Default is false, in which case we 
%                                     plot density as a function of 
%                                     eccentricity in visual degrees. Also
%                                     changes the reported density from per
%                                     deg^2 to per mm^2.
%
%    retinalMeridiansLegendsInsteadOfVisualSpaceMeridians - Logical.
%                                     Whether to label meridians in visual
%                                     space (default) or retinal space.
%
%
% Outputs:
%    hFig                           - The generated figure handle
%
% References:
%    Watson (2014). 'A formula for human RGC receptive field density as
%    a function of visual field location', JOV (2014), 14(7), 1-17.
%
% History:
%    11/8/19  NPC, ISETBIO Team     Wrote it.

    % Parse input
    p = inputParser;
    p.addParameter('eccentricityInMMInsteadOfDegs', false, @islogical);
    p.addParameter('whichEye','right',@(x)(ismember(x,{'left','right'})));
    p.addParameter('retinalMeridiansLegendsInsteadOfVisualSpaceMeridians', false, @islogical);
    p.parse(varargin{:});
     
    if (p.Results.eccentricityInMMInsteadOfDegs)
        eccentricity = obj.rhoDegsToMMs(obj.eccDegs);
    else
        eccentricity = obj.eccDegs;
    end
    whichEye = p.Results.whichEye;
    retinalMeridiansLegendsInsteadOfVisualSpaceMeridians = p.Results.retinalMeridiansLegendsInsteadOfVisualSpaceMeridians;
    
    % Reset figure
    if (isempty(hFig))
        hFig = figure();
    else
        figure(hFig); 
    end
    clf;
    legends = cell(numel(obj.enumeratedMeridianNames),1);
    
    % Label figure
    figureNumber = '1';
    set(hFig, 'Color', obj.figurePrefs.backgroundColor, 'Name', sprintf('Figure %s of %s', figureNumber, obj.paperTitleShort));
    
    % Loop over the four meridians
    for k = 1:numel(obj.enumeratedMeridianNames)
        
        % Get meridian name and color in Right Eye visual space
        meridianName = obj.enumeratedMeridianNames{k};
        theLegend = sprintf('%s (%s eye visual field)',meridianName, whichEye);
        theColor = obj.meridianColors(meridianName);
        
        % Correct meridian legend and color depending on whichEye and space of meridian (retinal or visual space)
        [theLegend, theColor] = obj.correctLegendAndColor(theLegend, theColor, meridianName, retinalMeridiansLegendsInsteadOfVisualSpaceMeridians, whichEye);
        
        % Accumulate legends
        legends{k} = theLegend;
        
        % Compute the data
        if (p.Results.eccentricityInMMInsteadOfDegs)
            eccUnits = 'mm';
            [~, coneRFDensity] = obj.coneRFSpacingAndDensity(eccentricity, meridianName, whichEye, eccUnits, 'Cones per mm2');
        else
            eccUnits = 'deg';
            [~, coneRFDensity] = obj.coneRFSpacingAndDensity(eccentricity, meridianName, whichEye, eccUnits, 'Cones per deg2');
        end
        
        % Plot the data
        plot(eccentricity, coneRFDensity, 'o-', ...
            'MarkerSize', 4, ...
            'MarkerFaceColor', min([1 1 1], theColor + 0.4), ...
            'Color', theColor, ...
            'LineWidth', obj.figurePrefs.markerLineWidth); 
        hold on;
    end
    
    % Add the maximal density line
    legends{numel(legends)+1} = '(0,0)';
    if (p.Results.eccentricityInMMInsteadOfDegs)
        plot([obj.eccDegs(2) obj.eccDegs(end)], obj.dc0/obj.alpha(0)*[1 1], 'k--', 'LineWidth', 1.5);
    else
        plot([obj.eccDegs(2) obj.eccDegs(end)], obj.dc0*[1 1], 'k--', 'LineWidth', 1.5);
    end
    
    % Show the legends
    legend(legends, 'Location', 'SouthWest');
    
    % Ticks and Lims
    if (p.Results.eccentricityInMMInsteadOfDegs)
        xlabel('eccentricity (retinal mm)', 'FontAngle', obj.figurePrefs.fontAngle);
        ylabel('density (cones/mm^2)', 'FontAngle', obj.figurePrefs.fontAngle);
        xLims = obj.rhoDegsToMMs([0.005 100]);
        yLims = [2000 250*1000];
        xTicks = [0.01 0.03 0.1 0.3 1 3 10];
        yTicks = [2000 5*1000 10*1000 20*1000 50*1000 100*1000 200*1000];
        yTicksLabels = {'2k', '5k', '10k', '20k', '50k', '100k', '200K'};
    else
        xlabel('eccentricity (degs)', 'FontAngle', obj.figurePrefs.fontAngle);
        ylabel('density (cones/deg^2)', 'FontAngle', obj.figurePrefs.fontAngle);
        xLims = [0.005 100];
        yLims = [100 20000];
        xTicks = [0.1 0.5 1 5 10 50 100];
        yTicks = [100 1000 10000];
        yTicksLabels = {'100', '1000', '10000'};
    end
    
    set(gca, 'XLim', xLims, 'YLim', yLims, ...
        'XScale', 'log', 'YScale', 'log', ...
        'XTick', xTicks, ...
        'YTick', yTicks, 'YTickLabel', yTicksLabels,...
        'FontSize', obj.figurePrefs.fontSize);
    
    % Grid & Title
    grid(gca, obj.figurePrefs.grid);
    title('Cone density as a function of eccentricity (Curcio et al, 1990)');
    drawnow;
end
