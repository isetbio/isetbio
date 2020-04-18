function unitTestFigure14(varargin)
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
    doIt(obj, eccDegs, eccUnits, densityUnits, meridianLabeling, 'mRGCToConesRatio', mfilename, plotlabOBJ);
end

function doIt(obj,eccentricities, eccUnits, densityUnits, meridianLabeling, figureName, theFileName, plotlabOBJ)

    exportFigure = false;
        
    hFig = figure(1); clf;
    theAxesGrid = plotlabOBJ.axesGrid(hFig, ...
            'leftMargin', 0.12, ...
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
        
        % cone RF spacing/density
        [coneRFSpacing, coneRFDensity] = obj.coneRFSpacingAndDensityAlongMeridian(eccentricities, ...
            rightEyeVisualFieldMeridianName, eccUnits, densityUnits);
       
        % mRGC spacing/density
        [mRGCRFSpacing, mRGCRFDensity] = obj.mRGCRFSpacingAndDensityAlongMeridian(eccentricities, ...
            rightEyeVisualFieldMeridianName, eccUnits, densityUnits);
        
        % ratios
        mRGCtoConeRatio = mRGCRFDensity./coneRFDensity;
        
        plot(eccentricities, mRGCtoConeRatio);
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
        xLabelString = sprintf('eccentricity (%s)', strrep(eccUnits, 'visual', ''));
    end
    
    yLims = [0.02 3];
    yTicks = [0.05 0.1 0.2 0.5 1 2];
    yTicksLabels = {'.05', '.10', '.20', '.50', '1.0', '2.0'};
    yLabelString = 'mRGC/cone ratio';
    
    
    % Labels and legends
    xlabel(xLabelString);
    ylabel(yLabelString);
    legend(theLegends, 'Location', 'SouthWest');
   
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
