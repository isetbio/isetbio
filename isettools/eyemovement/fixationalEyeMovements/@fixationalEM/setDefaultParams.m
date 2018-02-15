function setDefaultParams(obj)
    
    % Set the drift params from the Mergenthaler&Engbert, 2007 paper
    obj.setDriftParamsFromMergenthalerAndEngbert2007Paper();
    
    % Set microsaccade stats consistent with data in Martinez, 2009 
    obj.setMicroSaccadeStats();

    % How long to 'warm' up the drift model for
    obj.stabilizationSeconds = 2.0;
    
    % Microsaccade generation - related params
    obj.microSaccadeType = 'heatmap/fixation based';
    obj.heatMapWeight = 0.5;
    
    % Fixation map properties
    obj.fixationMapSpaceConstantArcMin = 15;

    % Heat map properties
    obj.heatMapWidthArcMin = 40;
    obj.heatMapSpatialSampleArcMin = 2;
    obj.heatMapTemporalSampleSeconds = 10/1000;
    obj.heatMapKernelSpaceConstantArcMin = 2.0;
    obj.heatMapKernelTimeConstantSeconds = 0.5;
    
    % Do not set the random seed
    obj.randomSeed = [];
    
    % Compute progress display
    obj.displayComputeProgress = false;
    obj.beVerbose = false;
end