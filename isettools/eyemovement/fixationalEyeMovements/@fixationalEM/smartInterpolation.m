function outMatrix = smartInterpolation(obj,inputTimeAxis, inputMatrix, outputTimeAxis)
    if (ndims(inputMatrix) > 3)
        error('input matrix must be 1D or 2D')
    elseif (ndims(inputMatrix)==3)
        if (length(inputTimeAxis) == size(inputMatrix,1))
            [T,N,M] = size(obj.heatMapTimeSeries);
            inputMatrix = reshape(inputMatrix,[T N*M]);
            outMatrix = reshape((interp1(inputTimeAxis, inputMatrix, outputTimeAxis)), [numel(outputTimeAxis) N M]);
        elseif (length(inputTimeAxis) == size(inputMatrix,3))
            [N,M,T] = size(obj.heatMapTimeSeries);
            inputMatrix = reshape(inputMatrix,[N*M T]);
            outMatrix = reshape((interp1(inputTimeAxis, inputMatrix', outputTimeAxis))', [numel(outputTimeAxis) N M]);
        else
            error('Time dimension should be either first or last')
        end
    else
        transposeOutMatrix = false;
        if (length(inputTimeAxis) == size(inputMatrix,2))
            inputMatrix = inputMatrix';
            transposeOutMatrix = true;
        end
        outMatrix = interp1(inputTimeAxis, inputMatrix, outputTimeAxis);
        if (transposeOutMatrix)
           outMatrix  = outMatrix'; 
        end
    end
end