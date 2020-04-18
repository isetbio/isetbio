function unitTestFigure11(varargin)
    % Parse input
    p = inputParser;
    p.addParameter('plotlabOBJ', [], @(x)(isempty(x) || isa(x, 'plotlab')));
    p.parse(varargin{:});
    plotlabOBJ = p.Results.plotlabOBJ;
    
    eccMinDegs = 0.0;
    eccMaxDegs = 10;
    eccSamplesNum = 200;
    eccDegs = linspace(eccMinDegs, eccMaxDegs, eccSamplesNum);
    eccUnits = 'deg';
    spacingUnits = 'deg';
    meridianLabeling = 'Watson'; %'retinal';   % choose from {'retinal', 'Watson'}
    
    obj = WatsonRGCModel();
    if (isempty(plotlabOBJ))
        plotlabOBJ  = obj.setUpPlotLab();
    end
    
    doIt(obj,eccDegs, eccUnits, spacingUnits , meridianLabeling, 'spacing', mfilename, plotlabOBJ);
end

function doIt(obj,eccentricities, eccUnits, spacingUnits, meridianLabeling, figureName, theFileName, plotlabOBJ)
    
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
        
        % One class On/OFF assuming equal numerosities
        OnOrOFF_mRGCRFSpacing = sqrt(2.0) * mRGCRFSpacing;
        OnOrOFF_mRGCRFDensity = 0.5 * mRGCRFDensity;
        
        if (strcmp(spacingUnits, 'deg'))
            % Convert to arcmin from deg
            OnOrOFF_mRGCRFSpacing = OnOrOFF_mRGCRFSpacing*60;
        else
            % Convert to microns from mm
            OnOrOFF_mRGCRFSpacing = OnOrOFF_mRGCRFSpacing*1000;
        end
        plot(theAxesGrid, eccentricities, OnOrOFF_mRGCRFSpacing);
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
        xLims = [0 10];
        xTicks = 0:2:10;
        xLabelString = sprintf('eccentricity (%s)', strrep(eccUnits, 'visual', ''));
    end
    
    if (strcmp(spacingUnits, 'mm'))
        yLims = [2000 250*1000];
        yTicks = [2000 5*1000 10*1000 20*1000 50*1000 100*1000 200*1000];
        yTicksLabels = {'2k', '5k', '10k', '20k', '50k', '100k', '200K'};
        yLabelString = sprintf('On or OFF spacing (microns)');
    else
        yLims = [0 6];
        yTicks = 0:1:6;
        yTicksLabels = 0:1:6;
        yLabelString = sprintf('spacing (arcmin)');
    end
    
    % Labels and legends
    xlabel(theAxesGrid, xLabelString);
    ylabel(theAxesGrid, yLabelString);
    legend(theAxesGrid, theLegends, 'Location', 'NorthWest');
   
    set(theAxesGrid, 'XLim', xLims, 'YLim', yLims, ...
        'XTick', xTicks, ...
        'YTick', yTicks, 'YTickLabel', yTicksLabels);
    
    % Export figure
    if (exportFigure)
        localDir = fileparts(which(theFileName));
        plotlabOBJ.exportFig(hFig, 'png', figureName, fullfile(localDir, 'exports'));
    end
end
