function unitTestFigure9()
    
    eccMinDegs = 0.1;
    eccMaxDegs = 80;
    eccSamplesNum = 100;
    eccDegs = logspace(log10(eccMinDegs), log10(eccMaxDegs), eccSamplesNum);
    eccUnits = 'deg';
    densityUnits = 'deg^2';
    meridianLabeling = 'Watson'; %'retinal';   % choose from {'retinal', 'Watson'}
    
    doIt(eccDegs, eccUnits, densityUnits, meridianLabeling);
end

function doIt(eccentricities, eccUnits, densityUnits, meridianLabeling)
    obj = WatsonRGCModel();
    meridianNames = obj.enumeratedMeridianNames;
    
    figure(1); clf; hold on;
    theLegends = cell(numel(meridianNames),1);
    
    % Loop over meridians
    for k = 1:numel(meridianNames)
        % Compute the data
        rightEyeVisualFieldMeridianName = meridianNames{k};
        [mRGCRFSpacing, mRGCRFDensity] = obj.mRGCRFSpacingAndDensityAlongMeridian(eccentricities, rightEyeVisualFieldMeridianName, eccUnits, densityUnits);
        
        plot(eccentricities, mRGCRFDensity, 'k-', ...
            'Color', obj.meridianColors(meridianNames{k}), ...
            'LineWidth', obj.figurePrefs.lineWidth);
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
        xLabelString = sprintf('eccentricity (%s)', strrep(eccUnits, 'visual', ''));
    end
    
    if (strcmp(densityUnits, 'retinal mm^2'))
        yLims = [2000 250*1000];
        yTicks = [2000 5*1000 10*1000 20*1000 50*1000 100*1000 200*1000];
        yTicksLabels = {'2k', '5k', '10k', '20k', '50k', '100k', '200K'};
        yLabelString = sprintf('density (mRGC RFs / %s)', densityUnits);
    else
        yLims = [1 40000];
        yTicks = [1 10 100 1000 10000];
        yTicksLabels = {'1' '10' '100', '1000', '10000'};
        yLabelString = sprintf('density (mRGC RFs / %s)', strrep(densityUnits, 'visual', ''));
    end
    
    % Labels and legends
    xlabel(xLabelString, 'FontAngle', obj.figurePrefs.fontAngle);
    ylabel(yLabelString, 'FontAngle', obj.figurePrefs.fontAngle);
    legend(theLegends);
   
    set(gca, 'XLim', xLims, 'YLim', yLims, ...
        'XScale', 'log', 'YScale', 'log', ...
        'XTick', xTicks, ...
        'YTick', yTicks, 'YTickLabel', yTicksLabels,...
        'FontSize', obj.figurePrefs.fontSize);
    
    % Finalize figure
    grid(gca, obj.figurePrefs.grid);
    
end
