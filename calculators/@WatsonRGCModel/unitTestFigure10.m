function unitTestFigure10(varargin)
    % Parse input
    p = inputParser;
    p.addParameter('plotlabOBJ', [], @(x)(isempty(x) || isa(x, 'plotlab')));
    p.parse(varargin{:});
    plotlabOBJ = p.Results.plotlabOBJ;
    
    eccMinDegs = 0.0;
    eccMaxDegs = 90;
    eccSamplesNum = 500;
    eccDegs = linspace(eccMinDegs, eccMaxDegs, eccSamplesNum);
    eccUnits = 'deg';
    spacingUnits = 'deg';
    meridianLabeling = 'Watson'; %'retinal';   % choose from {'retinal', 'Watson'}
    
    obj = WatsonRGCModel();
    if (isempty(plotlabOBJ))
        plotlabOBJ  = obj.setUpPlotLab();
    end
    
    doIt(obj, eccDegs, eccUnits, spacingUnits , meridianLabeling, 'spacing', mfilename, plotlabOBJ);
end

function doIt(obj, eccentricities, eccUnits, spacingUnits, meridianLabeling, figureName, theFileName, plotlabOBJ)
    exportFigure = false;
    
    hFig = figure(); clf;
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
        [mRGCRFSpacing, mRGCRFDensity] = obj.mRGCRFSpacingAndDensityAlongMeridian(eccentricities, rightEyeVisualFieldMeridianName, eccUnits, sprintf('%s^2',spacingUnits));
        plot(theAxesGrid, eccentricities, mRGCRFSpacing);
        if (strcmp(meridianLabeling, 'retinal'))
            theLegends{k} = rightEyeRetinalMeridianName;
        else
            theLegends{k} = rightEyeVisualFieldMeridianName;
        end
    end

     
    % Ticks and Lims
    if (strcmp(eccUnits, 'mm'))
        xLims = obj.rhoDegsToMMs([0.005 100]);
        xTicks = [0.01 0.03 0.1 0.3 1 3 10];
        xLabelString = sprintf('eccentricity (%s)', eccUnits);
    else
        xLims = [0 90];
        xTicks = 0:20:100;
        xLabelString = sprintf('eccentricity (%s)', strrep(eccUnits, 'visual', ''));
    end
    
    if (strcmp(spacingUnits, 'mm'))
        yLims = [2000 250*1000];
        yTicks = [2000 5*1000 10*1000 20*1000 50*1000 100*1000 200*1000];
        yTicksLabels = {'2k', '5k', '10k', '20k', '50k', '100k', '200K'};
        yLabelString = sprintf('spacing (microns)');
    else
        yLims = [0 1.02];
        yTicks = 0:0.2:1;
        yTicksLabels = 0:0.2:1;
        yLabelString = sprintf('spacing (deg)');
    end
    
    % Labels and legends
    xlabel(theAxesGrid, xLabelString);
    ylabel(theAxesGrid, yLabelString);
    legend(theAxesGrid, theLegends, 'Location', 'NorthWest')
   
    set(theAxesGrid,  'XLim', xLims, 'YLim', yLims, ...
        'XTick', xTicks, ...
        'YTick', yTicks, 'YTickLabel', yTicksLabels);
    
    % Export figure
    if (exportFigure)
        localDir = fileparts(which(theFileName));
        plotlabOBJ.exportFig(hFig, 'png', figureName, fullfile(localDir, 'exports'));
    end
    
end
