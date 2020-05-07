function fitModel(obj, varargin)
    p = inputParser;
    p.addParameter('dataset', 'medians', @(x)(ismember(x, {'medians', 'raw'})));
    p.parse(varargin{:});
    dataSet = p.Results.dataset;
    
    rng(1);
    
    if (strcmp(dataSet, 'raw'))
        fitRawData(obj);
    else
        fitMedianData(obj);
    end
    
    fitRetinalRawData(obj);
end

function fitRetinalRawData(obj)
    x = obj.centerData('size').retinalEccDegs;
    y = obj.centerData('size').retinalRadiusDegs;
    [obj.centerRetinalRadiusFunction, obj.centerRetinalRadiusParams, obj.centerRetinalRadiusParamsSE] = ...
        nonLinearFitData(x,y,[],[], obj.centerData('size').initialParams);
end


function fitRawData(obj)

    x = obj.centerData('size').eccDegs;
    y = obj.centerData('size').radiusDegs;
    [obj.centerRadiusFunction, obj.centerRadiusParams, obj.centerRadiusParamsSE] = ...
        nonLinearFitData(x,y,[],[], obj.centerData('size').initialParams);
    

    x = obj.surroundData('size').eccDegs;
    y = obj.surroundData('size').radiusDegs;
    [obj.surroundRadiusFunction, ...
        obj.surroundRadiusParams, obj.surroundRadiusParamsSE] = ...
        nonLinearFitData(x,y,[],[], obj.surroundData('size').initialParams);
    
    x = obj.centerData('sensitivity').radiusDegs;
    y = obj.centerData('sensitivity').peakSensitivity;
    [obj.centerPeakSensitivityFunction, ...
        obj.centerPeakSensitivityParams, obj.centerPeakSensitivityParamsSE] = ...
        nonLinearFitData(x,y,[],[], obj.centerData('sensitivity').initialParams);

    x = obj.surroundData('sensitivity').radiusDegs;
    y = obj.surroundData('sensitivity').peakSensitivity;
    [obj.surroundPeakSensitivityFunction,...
        obj.surroundPeakSensitivityParams, obj.surroundPeakSensitivityParamsSE] = ...
        nonLinearFitData(x,y,[],[], obj.surroundData('sensitivity').initialParams);
end


function fitMedianData(obj)
    % Fit the ecc - center radius data
    x = obj.centerData('size').eccDegsTable;
    yMedian = obj.centerData('size').radiusDegsMedianTable;
    yIQR = obj.centerData('size').radiusDegsIQRTable;
    ySamplesNum = obj.centerData('size').samplesTable;
    [obj.centerRadiusFunction, ...
        obj.centerRadiusParams, obj.centerRadiusParamsSE] = nonLinearFitData(...
        x,yMedian,yIQR,ySamplesNum, obj.centerData('size').initialParams);
    
    
    % Fit the ecc - surround radius data
    x = obj.surroundData('size').eccDegsTable;
    yMedian = obj.surroundData('size').radiusDegsMedianTable;
    yIQR = obj.surroundData('size').radiusDegsIQRTable;
    ySamplesNum = obj.surroundData('size').samplesTable;
    [obj.surroundRadiusFunction, ...
        obj.surroundRadiusParams, obj.surroundRadiusParamsSE] = nonLinearFitData(...
        x,yMedian,yIQR,ySamplesNum, obj.surroundData('size').initialParams);

    
    % Fit the radius - center sensitivity data
    x = obj.centerData('size').radiusDegsMedianTable;
    yMedian = obj.centerData('sensitivity').peakSensitivityMedianTable;
    yIQR = obj.centerData('sensitivity').peakSensitivityIQRTable;
    ySamplesNum = obj.centerData('sensitivity').samplesTable;
    [obj.centerPeakSensitivityFunction, ...
        obj.centerPeakSensitivityParams, obj.centerPeakSensitivityParamsSE] = nonLinearFitData(...
        x,yMedian,yIQR,ySamplesNum, obj.centerData('sensitivity').initialParams);
    
    
    % Fit the radius - surround sensitivity data
    x = obj.surroundData('size').radiusDegsMedianTable;
    yMedian = obj.surroundData('sensitivity').peakSensitivityMedianTable;
    yIQR = obj.surroundData('sensitivity').peakSensitivityIQRTable;
    ySamplesNum = obj.surroundData('sensitivity').samplesTable;
    [obj.surroundPeakSensitivityFunction, ...
        obj.surroundPeakSensitivityParams, obj.surroundPeakSensitivityParamsSE] = nonLinearFitData(...
        x,yMedian,yIQR,ySamplesNum, obj.surroundData('sensitivity').initialParams); 
end


function [powerFunction, fittedParams, fittedParamsSE] = nonLinearFitData(x,y, yIQR,ySamplesNum,initialParams)
    powerFunction = @(p,x)(p(1)*x.^p(2));  % Objective Function
    opts.RobustWgtFun = 'talwar';
    if (~isempty(yIQR))
        xx = [];
        yy = [];
        for k = 1:numel(y)
            ySamples = normrnd(y(k), yIQR(k)/1.35, [1 ySamplesNum(k)]);
            yy = cat(2, yy, ySamples);
            xx = cat(2, xx, repmat(x(k), [1 ySamplesNum(k)]));
        end
        y = yy;
        x = xx;
    end
    [fittedParams,~,~,varCovarianceMatrix,~] = nlinfit(x,y,powerFunction,initialParams,opts);
    % standard error of the mean
    fittedParamsSE = sqrt(diag(varCovarianceMatrix));
    fittedParamsSE = fittedParamsSE';
   
    % make it standard deviation
    %fittedParamsSE = fittedParamsSE * sqrt(mean(ySamplesNum));
end
