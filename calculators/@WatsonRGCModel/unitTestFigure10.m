function unitTestFigure10()
    
    eccMinDegs = 0.0;
    eccMaxDegs = 90;
    eccSamplesNum = 500;
    eccDegs = linspace(eccMinDegs, eccMaxDegs, eccSamplesNum);
    eccUnits = 'deg';
    spacingUnits = 'deg';
    meridianLabeling = 'Watson'; %'retinal';   % choose from {'retinal', 'Watson'}
    
    doIt(eccDegs, eccUnits, spacingUnits , meridianLabeling);
end

function doIt(eccentricities, eccUnits, spacingUnits, meridianLabeling)
    obj = WatsonRGCModel();
    meridianNames = obj.enumeratedMeridianNames;
    
    figure(1); clf; hold on;
    theLegends = cell(numel(meridianNames),1);
    
    % Loop over meridians
    for k = 1:numel(meridianNames)
        % Compute the data
        rightEyeVisualFieldMeridianName = meridianNames{k};
        [mRGCRFSpacing, mRGCRFDensity] = obj.mRGCRFSpacingAndDensityAlongMeridian(eccentricities, rightEyeVisualFieldMeridianName, eccUnits, sprintf('%s^2',spacingUnits));
        if (strcmp(spacingUnits, 'deg'))
            % Convert to arcmin from deg
            %mRGCRFSpacing = mRGCRFSpacing*60;
        else
            % Convert to microns from mm
            %mRGCRFSpacing = mRGCRFSpacing*1000;
        end
        plot(eccentricities, mRGCRFSpacing, 'k-', ...
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
    xlabel(xLabelString, 'FontAngle', obj.figurePrefs.fontAngle);
    ylabel(yLabelString, 'FontAngle', obj.figurePrefs.fontAngle);
    legend(theLegends);
   
    set(gca, 'XLim', xLims, 'YLim', yLims, ...
        'XScale', 'linear', 'YScale', 'linear', ...
        'XTick', xTicks, ...
        'YTick', yTicks, 'YTickLabel', yTicksLabels,...
        'FontSize', obj.figurePrefs.fontSize);
    
    % Finalize figure
    grid(gca, obj.figurePrefs.grid);
    
end
