function plotDigitizedData(obj, varargin)
    p = inputParser;
    p.addParameter('generatePlots', true, @islogical);
    p.parse(varargin{:});
    generatePlots = p.Results.generatePlots;
    hFig = figure(1); clf;
    
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'rowsNum', 2, ...
            'colsNum', 2, ...
            'leftMargin', 0.08, ...
            'rightMargin', 0.01, ...
            'widthMargin', 0.1, ...
            'heightMargin', 0.06, ...
            'bottomMargin', 0.06, ...
            'topMargin', 0.04);
      
    eccDegs = logspace(log10(0.1), log10(180), 100);
    
    centerRadiusModel.eccDegs = eccDegs;
    centerRadiusModel.function = obj.centerRadiusFunction;
    centerRadiusModel.params = obj.centerRadiusParams;
    
    surroundRadiusModel.eccDegs = eccDegs;
    surroundRadiusModel.function = obj.surroundRadiusFunction;
    surroundRadiusModel.params = obj.surroundRadiusParams;
    
    centerPeakSensitivityModel.radiusDegs = centerRadiusModel.function(centerRadiusModel.params, eccDegs);
    centerPeakSensitivityModel.function = obj.centerPeakSensitivityFunction;
    centerPeakSensitivityModel.params = obj.centerPeakSensitivityParams;
    
    surroundPeakSensitivityModel.radiusDegs = surroundRadiusModel.function(surroundRadiusModel.params, eccDegs);
    surroundPeakSensitivityModel.function = obj.surroundPeakSensitivityFunction;
    surroundPeakSensitivityModel.params = obj.surroundPeakSensitivityParams;
    
    pointSize = 144;
    CronerKaplanRGCModel.plotRadii(theAxesGrid{1,1}, obj.centerData, centerRadiusModel, pointSize, [1 0 0], true, 'center');
    CronerKaplanRGCModel.plotRadii(theAxesGrid{1,2}, obj.surroundData, surroundRadiusModel, pointSize, [0 0 1], false, 'surround');
    
    CronerKaplanRGCModel.plotSensitivities(theAxesGrid{2,1}, obj.centerData, centerPeakSensitivityModel, pointSize, [1 0 0], true, '');
    CronerKaplanRGCModel.plotSensitivities(theAxesGrid{2,2}, obj.surroundData, surroundPeakSensitivityModel, pointSize, [0 0 1], false, '');

    if (generatePlots)
        obj.plotlabOBJ.exportFig(hFig, 'pdf', 'CronnerKaplanParams', pwd());
    end
    
    % Plot the retinal radius
    hFig = figure(11); clf
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'rowsNum', 2, ...
            'colsNum', 2, ...
            'leftMargin', 0.08, ...
            'rightMargin', 0.01, ...
            'widthMargin', 0.1, ...
            'heightMargin', 0.06, ...
            'bottomMargin', 0.06, ...
            'topMargin', 0.04);
    

    centerRetinalRadiusModel.eccDegs = eccDegs;
    centerRetinalRadiusModel.function = obj.centerRetinalRadiusFunction;
    centerRetinalRadiusModel.params = obj.centerRetinalRadiusParams;
    CronerKaplanRGCModel.plotRadii(theAxesGrid{1,1}, obj.centerData, centerRetinalRadiusModel, pointSize, [1 0 0], true, 'retinal center');

    eccMM = 0.3:0.1:20;
    eccDegsFromMM = @(x)(0.1 + 4.21 * x + 0.038 * x.^2);
    midgetDendriticFieldDiameterMM =  8.64 * eccMM.^1.04/1000;
    midgetDendriticFieldDiameterDegs = eccDegsFromMM(midgetDendriticFieldDiameterMM);
    
    
    plotRadiusRatio(theAxesGrid{1,2}, centerRadiusModel, eccDegsFromMM(eccMM), midgetDendriticFieldDiameterDegs, [1 0 0]);
    pause
end

function plotRadiusRatio(theAxes, visualModel, retinalModelEccDegs, retinalModelDendriticDiameter,  color)

    hold(theAxes, 'on');
    
    visualRadius  = visualModel.function(visualModel.params, visualModel.eccDegs);
    size(retinalModelEccDegs)
    size(retinalModelDendriticDiameter)
    line(theAxes, retinalModelEccDegs, retinalModelDendriticDiameter, 'Color', [0 0 0]);
    line(theAxes, visualModel.eccDegs, visualRadius, 'Color', color);
    
    set(theAxes, 'XLim', [0 100], 'YLim', [0 0.3]);
    set(theAxes, 'XTick', 0:20:100, 'YTick', 0:0.1:0.3);
    xlabel(theAxes, 'eccentricity (degs)');
    ylabel(theAxes, 'radius');
    set(theAxes, 'XScale', 'linear', 'YScale', 'linear');
    grid(theAxes, 'on');
    box(theAxes, 'off');
    
end
