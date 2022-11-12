function  c = weightedMean(data, weights)
    
    dataPointsNum = size(data,1);
    dimensionsNum = size(data,2);
    assert(dataPointsNum == numel(weights), 'Inconsistent data');
    
    c = zeros(1, dimensionsNum);
    w = weights(:);
    ss = sum(w);
    
    for iDim = 1:dimensionsNum
        c(:, iDim) = dot(data(:, iDim), w,1)/ss;
    end
end