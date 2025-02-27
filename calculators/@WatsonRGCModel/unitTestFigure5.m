function plotlabOBJ = unitTestFigure5(varargin)
% Generate Figure 5 of Watson (2014) which plots total RGC RF density as a function
% of eccentricity.

    % Parse input
    p = inputParser;
    p.addParameter('plotlabOBJ', [], @(x)(isempty(x) || isa(x, 'plotlab')));
    p.parse(varargin{:});
    plotlabOBJ = p.Results.plotlabOBJ;
    
    eccMinDegs = 0.1;
    eccMaxDegs = 100;
    eccSamplesNum = 100;
    eccDegs = logspace(log10(eccMinDegs), log10(eccMaxDegs), eccSamplesNum);
    eccUnits = 'deg';
    densityUnits = 'deg^2';
    meridianLabeling = 'Watson'; %'retinal';   % choose from {'retinal', 'Watson'}
    
    obj = WatsonRGCModel();
    if (isempty(plotlabOBJ))
        plotlabOBJ  = obj.setUpPlotLab();
    end
    
    doIt(obj, eccDegs, eccUnits, densityUnits, meridianLabeling, 'TotalRGCDensity', mfilename, plotlabOBJ);
end

function doIt(obj, eccentricities, eccUnits, densityUnits, meridianLabeling, figureName, theFileName, plotlabOBJ)

    exportFigure = false;
        
    hFig = figure(5); clf;
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
        % Compute the data
        rightEyeVisualFieldMeridianName = meridianNames{k};
        totalRGCRFDensity = obj.totalRGCRFDensityAlongMeridian(eccentricities, rightEyeVisualFieldMeridianName, eccUnits, densityUnits);
        
        plot(eccentricities, totalRGCRFDensity);
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
        xTicks = [0.1 0.3 1 3 10 30 100];
        totalConesAtFovea = obj.dc0;
        midgetFractionAtFovea = 1/obj.f0;
        totalRGCsAtFovea = totalConesAtFovea * 2 / midgetFractionAtFovea;
        plot(theAxesGrid, xLims(1), totalRGCsAtFovea, 'ko');
        xLabelString = sprintf('eccentricity (%s)', strrep(eccUnits, 'visual', ''));
    end
    
    if (strcmp(densityUnits, 'retinal mm^2'))
        yLims = [2000 250*1000];
        yTicks = [2000 5*1000 10*1000 20*1000 50*1000 100*1000 200*1000];
        yTicksLabels = {'2k', '5k', '10k', '20k', '50k', '100k', '200K'};
        yLabelString = sprintf('density (mRGC RFs / %s)', densityUnits);
    else
        yLims = [4 40000];
        yTicks = [10 100 1000 10000];
        yTicksLabels = {'10' '100', '1000', '10000'};
        yLabelString = sprintf('density (total RGC RFs / %s)', strrep(densityUnits, 'visual', ''));
    end
    
    % Labels and legends
    xlabel(xLabelString);
    ylabel(yLabelString);
    legend(theLegends);
   
    set(gca, 'XLim', xLims, 'YLim', yLims, ...
        'XScale', 'log', 'YScale', 'log', ...
        'XTick', xTicks, ...
        'YTick', yTicks, 'YTickLabel', yTicksLabels);
    
    % Export figure
    if (exportFigure)
        localDir = fileparts(which(theFileName));
        plotlabOBJ.exportFig(hFig, 'png', figureName, fullfile(localDir, 'exports'));
    end
end
