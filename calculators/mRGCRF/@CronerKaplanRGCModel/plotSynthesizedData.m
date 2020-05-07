function [hFig1, hFig2, hFig3, hFig4] = plotSynthesizedData(obj)

    % Find appropriate scatter point size
    dataSampleNum = numel(obj.synthesizedData.eccDegs);
    if (dataSampleNum < 200)
        pointSize = 144;
    else
        r = 2500/dataSampleNum;
        if (r < 1)
            r = 1;
        end
        pointSize = r^2;
    end
    
    eccDegsDisplayed = prctile(obj.synthesizedData.eccDegs, [5 10 20 40]);
    
    [hFig3, hFig4] = plotRFinstances(obj, eccDegsDisplayed);
    hFig2 = plotSynthesizedRatios(obj, pointSize);
    hFig1 = plotSynthesizedRadiiAndSensitivities(obj, pointSize);
    p = get(hFig1, 'Position');
    pRight = p;
    pRight(1) = 14.3;
    set(hFig1, 'Position', pRight);
end


function hFig = plotSynthesizedRadiiAndSensitivities(obj, pointSize)
    % Retrieve the data
    eccDegs = obj.synthesizedData.eccDegs;
    centerRadii = obj.synthesizedData.centerRadii;
    surroundRadii = obj.synthesizedData.surroundRadii; 
    centerPeakSensitivities = obj.synthesizedData.centerPeakSensitivities;
    surroundPeakSensitivities = obj.synthesizedData.surroundPeakSensitivities;
     
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
    
    centerData = containers.Map();
    centerData('size') = struct(...
        'eccDegs', eccDegs, ...
        'radiusDegs', centerRadii);
    centerData('sensitivity') = struct(...
        'radiusDegs', centerRadii, ...
        'peakSensitivity',  centerPeakSensitivities);
    
    surroundData = containers.Map();
    surroundData('size') = struct(...
        'eccDegs', eccDegs, ...
        'radiusDegs', surroundRadii);
    surroundData('sensitivity') = struct(...
        'radiusDegs', surroundRadii, ...
        'peakSensitivity',  surroundPeakSensitivities);
    
    
    hFig = figure(3); clf;
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'rowsNum', 2, ...
            'colsNum', 2, ...
            'leftMargin', 0.08, ...
            'rightMargin', 0.01, ...
            'widthMargin', 0.1, ...
            'heightMargin', 0.06, ...
            'bottomMargin', 0.06, ...
            'topMargin', 0.04);
      
    
    CronerKaplanRGCModel.plotRadii(theAxesGrid{1,1}, centerData, centerRadiusModel, pointSize,[1 0 0], true, 'center');
    CronerKaplanRGCModel.plotRadii(theAxesGrid{1,2}, surroundData, surroundRadiusModel, pointSize,[0 0 1], false, 'surround');
    
    CronerKaplanRGCModel.plotSensitivities(theAxesGrid{2,1}, centerData, centerPeakSensitivityModel, pointSize,[1 0 0], true, '');
    CronerKaplanRGCModel.plotSensitivities(theAxesGrid{2,2}, surroundData, surroundPeakSensitivityModel, pointSize,[0 0 1], false, '');

end

function hFig = plotSynthesizedRatios(obj, pointSize)
    
    % Retrieve the data
    eccDegs = obj.synthesizedData.eccDegs;
    centerRadii = obj.synthesizedData.centerRadii;
    surroundRadii = obj.synthesizedData.surroundRadii; 
    centerPeakSensitivities =obj.synthesizedData.centerPeakSensitivities;
    surroundPeakSensitivities = obj.synthesizedData.surroundPeakSensitivities;
    
    centerIntegratedSensitivities = pi * centerRadii.^2 .* centerPeakSensitivities;
    surroundIntegratedSensitivities = pi * surroundRadii.^2 .* surroundPeakSensitivities;
    
    hFig = figure(2); clf;
    
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'rowsNum', 2, ...
            'colsNum', 2, ...
            'leftMargin', 0.08, ...
            'rightMargin', 0.01, ...
            'widthMargin', 0.1, ...
            'heightMargin', 0.06, ...
            'bottomMargin', 0.06, ...
            'topMargin', 0.04);
        
    % Ratio of center/surround radius as a function of eccentricity
    theAxes = theAxesGrid{1,1};
    ratios = centerRadii./surroundRadii;
    scatter(theAxes, eccDegs, ratios, pointSize);
    % Fit a line to the ratios
    initialParams = [0.1 1];
    [xFit,yFit] = fitLinearFunction(eccDegs,ratios, initialParams);
    % plot the fitted line
    hold(theAxes, 'on');
    plot(theAxes,xFit, yFit, 'k-');
    set(theAxes, 'XLim', [0 40], 'XTick', 0:10:40, 'YLim', [0 1], 'YTick', 0:0.1:1);
    xlabel(theAxes, 'eccentricity (degs)');
    ylabel(theAxes, 'radius (center:surround)');
    grid(theAxes, 'on');
    box(theAxes, 'off');
    
    % Ratio of surround/center peak sensitivity as a function of eccentricity
    theAxes = theAxesGrid{1,2};
    ratios = surroundPeakSensitivities./centerPeakSensitivities;
    scatter(theAxes, eccDegs, ratios, pointSize);
%     % Fit a line to the ratios
%     initialParams = [0.01 1];
%     [xFit,yFit] = fitLinearFunction(eccDegs,ratios, initialParams);
%     % plot the fitted line
%     hold(theAxes, 'on');
%     plot(theAxes,xFit, yFit, 'k-');
    set(theAxes, 'XLim', [0.1 100], 'XTick', [0.1 1 10 100]);
    set(theAxes, 'XScale', 'log', 'YScale', 'log', 'YLim', [0.0001 1], 'YTick', [0.0001 0.001 0.01 0.1 1]);
    xlabel(theAxes, 'eccentricity (degs)');
    ylabel(theAxes, 'peak sensitivity (surround:center)');
    grid(theAxes, 'on');
    box(theAxes, 'off');
    
    % Ratio of surround/center integrated sensitivity as a function of eccentricity
    theAxes = theAxesGrid{2,1};
    ratios = surroundIntegratedSensitivities./centerIntegratedSensitivities;
    scatter(theAxes, eccDegs, ratios, pointSize);
    % Fit a line to the ratios
    initialParams = [0.05 1];
    [xFit,yFit] = fitLinearFunction(eccDegs,ratios, initialParams);
    % plot the fitted line
    hold(theAxes, 'on');
    plot(theAxes,xFit, yFit, 'k-');
    set(theAxes, 'XLim', [0 40], 'XTick', 0:10:40);
    set(theAxes, 'YLim', [0 1.0], 'YTick', 0:0.2:5);
    xlabel(theAxes, 'eccentricity (degs)');
    ylabel(theAxes, 'integrated sensitivity (surround:center)');
    grid(theAxes, 'on');
    box(theAxes, 'off');
    
    % surround/center integrated sensitivity 
    theAxes = theAxesGrid{2,2};
    scatter(theAxes, centerIntegratedSensitivities, surroundIntegratedSensitivities, pointSize);
%     % Fit a line 
%     initialParams = [0.0 1];
%     [xFit,yFit] = fitLinearFunction(centerIntegratedSensitivities, surroundIntegratedSensitivities, initialParams);
%     % plot the fitted line
%     hold(theAxes, 'on');
%     plot(theAxes,xFit, yFit, 'k-');
    set(theAxes, 'XLim', [0 4], 'XTick', 0:1:12);
    set(theAxes, 'YLim', [0 4], 'YTick', 0:1:12);
    xlabel(theAxes, 'integrated sensitivity (center)');
    ylabel(theAxes, 'integrated sensitivity (surround)');
    grid(theAxes, 'on');
    box(theAxes, 'off');
end

function [xFit, yFit] = fitLinearFunction(x,y, initialParams)
    linearFunction = @(p,x)(p(1) + p(2)*x);  % Objective Function
    costFunction = @(p) norm(y - linearFunction(p,x));     % Residual Norm Cost Function
    fitParams = fminsearch(costFunction, initialParams);
    xFit = linspace(min(x), max(x), 100);
    yFit = fitParams(1) + fitParams(2)*xFit;
end
