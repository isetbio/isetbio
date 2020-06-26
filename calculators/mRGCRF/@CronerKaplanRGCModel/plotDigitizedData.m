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
end

