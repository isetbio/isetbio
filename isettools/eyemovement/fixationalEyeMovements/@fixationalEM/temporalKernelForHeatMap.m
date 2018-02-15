% Method to compute the temporal kernel for the heat map
function [kernel, kernelSupport] = temporalKernelForHeatMap(obj, gridNodesNum)
    kernelLength = round(obj.heatMapKernelTimeConstantSeconds/(obj.heatMapUpdateIntervalStepsNum*obj.timeStepDurationSeconds));
    kernelSupport = (0:(kernelLength-1))*obj.heatMapUpdateIntervalStepsNum*obj.timeStepDurationSeconds;
    if (obj.heatMapKernelTimeConstantSeconds < obj.timeStepDurationSeconds)
        kernelTimeConstantSeconds = obj.timeStepDurationSeconds;
    else
        kernelTimeConstantSeconds = obj.heatMapKernelTimeConstantSeconds;
    end
    kernel = exp(-(kernelSupport/kernelTimeConstantSeconds));
    
    % subsample according to obj.heatMapTemporalSampleSeconds
    indices = 1:round(kernelSupport(end)/obj.heatMapTemporalSampleSeconds);
    kernel = kernel(indices);
    kernelSupport = kernelSupport(indices);
    % zero kernel below 1%
    kernel = kernel(kernel>=0.01);
    kernel = (kernel-min(kernel))/(max(kernel)-min(kernel));
    % flip kernel because instead of convolution we are doing a dot product
    kernel = fliplr(kernel(indices));
    kernel = kernel';
    kernelSupport = kernelSupport';
    kernel = repmat(kernel, [1  gridNodesNum gridNodesNum]);
end