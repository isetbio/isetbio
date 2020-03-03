function unitTestFigure14()
    
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
        
        % cone RF spacing/density
        [coneRFSpacing, coneRFDensity] = obj.coneRFSpacingAndDensityAlongMeridian(eccentricities, ...
            rightEyeVisualFieldMeridianName, eccUnits, densityUnits);
       
        % mRGC spacing/density
        [mRGCRFSpacing, mRGCRFDensity] = obj.mRGCRFSpacingAndDensityAlongMeridian(eccentricities, ...
            rightEyeVisualFieldMeridianName, eccUnits, densityUnits);
        
        % ratios
        mRGCtoConeRatio = mRGCRFDensity./coneRFDensity;
        
        plot(eccentricities, mRGCtoConeRatio, 'k-', ...
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
        xLims = [0.2 80];
        xTicks = [0.5 1 5 10 50];
        xLabelString = sprintf('eccentricity (%s)', strrep(eccUnits, 'visual', ''));
    end
    
    yLims = [0.02 2.2];
    yTicks = [0.05 0.1 0.2 0.5 1 2];
    yTicksLabels = {'.05', '.10', '.20', '.50', '1.0', '2.0'};
    yLabelString = 'mRGC/cone ratio';
    
    
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

function old()
    meridianNames = obj.enumeratedMeridianNames;
    
    figure(1); clf; hold on;
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
        plot(eccentricities, OnOrOFF_mRGCRFSpacing, 'k-', ...
            'Color', obj.meridianColors(meridianNames{k}), ...
            'LineWidth', obj.figurePrefs.lineWidth);
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
    xlabel(xLabelString, 'FontAngle', obj.figurePrefs.fontAngle);
    ylabel(yLabelString, 'FontAngle', obj.figurePrefs.fontAngle);
    legend(theLegends, 'Location', 'NorthWest');
   
    set(gca, 'XLim', xLims, 'YLim', yLims, ...
        'XScale', 'linear', 'YScale', 'linear', ...
        'XTick', xTicks, ...
        'YTick', yTicks, 'YTickLabel', yTicksLabels,...
        'FontSize', obj.figurePrefs.fontSize);
    
    % Finalize figure
    grid(gca, obj.figurePrefs.grid);
    
end
