function plotlabOBJ = unitTestFigure1(varargin)
% Generate Figure 1 of Watson (2014) which plots cone density as a function
% of eccentricity using Curcio's 1990 data.

    % Parse input
    p = inputParser;
    p.addParameter('plotlabOBJ', [], @(x)(isempty(x) || isa(x, 'plotlab')));
    p.parse(varargin{:});
    plotlabOBJ = p.Results.plotlabOBJ;
    
    eccMinDegs = 0.2;
    eccMaxDegs = 100;
    eccSamplesNum = 100;
    eccDegs = logspace(log10(eccMinDegs), log10(eccMaxDegs), eccSamplesNum);
    eccUnits = 'deg';
    densityUnits = 'deg^2';
    meridianLabeling = 'Watson'; %'retinal';   % choose from 'retinal', 'Watson'

    obj = WatsonRGCModel();
    if (isempty(plotlabOBJ))
        plotlabOBJ  = obj.setUpPlotLab();
    end
    
    doIt(obj, eccDegs, eccUnits, densityUnits, meridianLabeling, 'coneDensity', mfilename, plotlabOBJ);
end

function doIt(obj, eccentricities, eccUnits, densityUnits, meridianLabeling, figureName, theFileName, plotlabOBJ)
    
    exportFigure = false;
    
    hFig = figure(1); clf;
    theAxesGrid = plotlabOBJ.axesGrid(hFig, ...
            'leftMargin', 0.16, ...
            'bottomMargin', 0.18, ...
            'rightMargin', 0.04, ...
            'topMargin', 0.05);
    theAxesGrid = theAxesGrid{1,1};
    hold(theAxesGrid, 'on');
    
    meridianNames = obj.enumeratedMeridianNames;
    theLegends = cell(numel(meridianNames),1);
    
    % Loop over meridians
    for k = 1:numel(meridianNames)
        rightEyeVisualFieldMeridianName = meridianNames{k};
        [coneRFSpacing, coneRFDensity, rightEyeRetinalMeridianName] = ...
            obj.coneRFSpacingAndDensityAlongMeridian(eccentricities, rightEyeVisualFieldMeridianName, eccUnits, densityUnits);
        plot(theAxesGrid, eccentricities, coneRFDensity);
        if (strcmp(meridianLabeling, 'retinal'))
            theLegends{k} = rightEyeRetinalMeridianName;
        else
            theLegends{k} = rightEyeVisualFieldMeridianName;
        end
    end

     
    % Ticks and Lims
    if (strcmp(eccUnits, 'retinal mm'))
        xLims = obj.rhoDegsToMMs([0.005 100]);
        xTicks = [0.01 0.03 0.1 0.3 1 3 10];
        xLabelString = sprintf('eccentricity (%s)', eccUnits);
    else
        xLims = [0.05 100];
        xTicks = [0.1 0.5 1 5 10 50 100];
        plot(theAxesGrid, xLims(1), obj.dc0, 'ko');
        xLabelString = sprintf('eccentricity (%s)', strrep(eccUnits, 'visual', ''));
    end
    
    if (strcmp(densityUnits, 'retinal mm^2'))
        yLims = [2000 250*1000];
        yTicks = [2000 5*1000 10*1000 20*1000 50*1000 100*1000 200*1000];
        yTicksLabels = {'2k', '5k', '10k', '20k', '50k', '100k', '200K'};
        yLabelString = sprintf('density (cones / %s)', densityUnits);
    else
        yLims = [100 20000];
        yTicks = [1 10 100 1000 10000];
        yTicksLabels = {'1' '10' '100', '1000', '10000'};
        yLabelString = sprintf('density (cones / %s)', strrep(densityUnits, 'visual', ''));
    end
    
    % Labels and legends
    xlabel(theAxesGrid, xLabelString);
    ylabel(theAxesGrid, yLabelString);
    legend(theAxesGrid,theLegends);
   
    set(theAxesGrid, 'XLim', xLims, 'YLim', yLims, ...
        'XScale', 'log', 'YScale', 'log', ...
        'XTick', xTicks, 'MinorGridLineStyle', '-', 'MinorGridAlpha', 0.1, ...
        'YTick', yTicks, 'YTickLabel', yTicksLabels);
    
    % Export figure
    if (exportFigure)
        localDir = fileparts(which(theFileName));
        plotlabOBJ.exportFig(hFig, 'png', figureName, fullfile(localDir, 'exports'));
    end
end
