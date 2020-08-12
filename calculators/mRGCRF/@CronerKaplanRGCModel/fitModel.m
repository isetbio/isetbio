function fitModel(obj, varargin)
    p = inputParser;
    p.addParameter('dataset', 'medians', @(x)(ismember(x, {'medians', 'raw', 'paperFormulas'})));
    p.parse(varargin{:});
    dataSet = p.Results.dataset;
    
    rng(1);
    
    switch (dataSet)
        case 'raw'
            fitRawData(obj);
        case 'medians'
            fitMedianData(obj);
        case 'paperFormulas'
            usePaperFits(obj)
        otherwise
            error('Unknown dataSet: ''%s''.', dataSet)
    end

end


function fitRawData(obj)

    x = obj.centerData('size').eccDegs;
    y = obj.centerData('size').radiusDegs;
    obj.centerRadiusThreshold = min(obj.centerData('size').radiusDegs);
    [obj.centerRadiusFunction, obj.centerRadiusParams, obj.centerRadiusParamsSE] = ...
        nonLinearFitDataWithThreshold(x,y, obj.centerRadiusThreshold, [],[], obj.centerData('size').initialParams);

    x = obj.surroundData('size').eccDegs;
    y = obj.surroundData('size').radiusDegs;
    obj.surroundRadiusThreshold = min(obj.surroundData('size').radiusDegs);
    [obj.surroundRadiusFunction, ...
        obj.surroundRadiusParams, obj.surroundRadiusParamsSE] = ...
        nonLinearFitDataWithThreshold(x,y,obj.surroundRadiusThreshold,[],[], obj.surroundData('size').initialParams);
    
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

function usePaperFits(obj)

    % Fit the ecc - center radius data
    x = obj.centerData('size').eccDegsTable;
    yMedian = obj.centerData('size').radiusDegsMedianTable;
    yIQR = obj.centerData('size').radiusDegsIQRTable;
    ySamplesNum = obj.centerData('size').samplesTable;
    [obj.centerRadiusFunction, ...
        obj.centerRadiusParams, obj.centerRadiusParamsSE] = nonLinearFitData(...
        x,yMedian,yIQR,ySamplesNum, obj.centerData('size').initialParams);
    
    % The ecc - surround radius equation from the paper (Figure 4 caption)
    obj.surroundRadiusFunction = @(p,x)(p(1)*x.^p(2));
    obj.surroundRadiusParams = [0.203 0.472];
    obj.surroundRadiusParamsSE = [0 0];
    
    % The radius - center sensitivity equation from the paper (Figure 5 caption)
    obj.centerPeakSensitivityFunction = @(p,x)(p(1)*x.^p(2));
    obj.centerPeakSensitivityParams = [0.391 -1.850];
    obj.centerPeakSensitivityParamsSE = [0 0];
    
    % The radius - surround sensitivity equation from the paper (Figure 5 caption)
    obj.surroundPeakSensitivityFunction = @(p,x)(p(1)*x.^p(2));
    obj.surroundPeakSensitivityParams = [0.128 -2.147];
    obj.surroundPeakSensitivityParamsSE = [0 0];
end

function fitMedianData(obj)
    % Fit the ecc - center radius data
    x = obj.centerData('size').eccDegsTable;
    yMedian = obj.centerData('size').radiusDegsMedianTable;
    yIQR = obj.centerData('size').radiusDegsIQRTable;
    ySamplesNum = obj.centerData('size').samplesTable;
    obj.centerRadiusThreshold = min(obj.centerData('size').radiusDegs);
    [obj.centerRadiusFunction, ...
        obj.centerRadiusParams, obj.centerRadiusParamsSE] = nonLinearFitDataWithThreshold(...
        x, yMedian, obj.centerRadiusThreshold, yIQR,ySamplesNum, obj.centerData('size').initialParams);
    
    % Fit the ecc - surround radius data
    x = obj.surroundData('size').eccDegsTable;
    yMedian = obj.surroundData('size').radiusDegsMedianTable;
    yIQR = obj.surroundData('size').radiusDegsIQRTable;
    ySamplesNum = obj.surroundData('size').samplesTable;
    obj.surroundRadiusThreshold = min(obj.surroundData('size').radiusDegs);
    [obj.surroundRadiusFunction, ...
        obj.surroundRadiusParams, obj.surroundRadiusParamsSE] = nonLinearFitDataWithThreshold(...
        x,yMedian, obj.surroundRadiusThreshold, yIQR,ySamplesNum, obj.surroundData('size').initialParams);

    
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

function [powerFunction,  fittedParams, fittedParamsSE] = nonLinearFitDataWithThreshold(x,y, threshold, yIQR,ySamplesNum,initialParams)
    powerFunction = @(p, x)( max(threshold*ones(size(x)), p(1)*(abs(x)).^p(2)) );  % Objective Function
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
    [fittedParams,~,~,varCovarianceMatrix,~] = nlinfit(x,y,powerFunction,initialParams, opts);
    % standard error of the mean
    fittedParamsSE = sqrt(diag(varCovarianceMatrix));
    fittedParamsSE = fittedParamsSE';
   
    % make it standard deviation
    %fittedParamsSE = fittedParamsSE * sqrt(mean(ySamplesNum));
end

