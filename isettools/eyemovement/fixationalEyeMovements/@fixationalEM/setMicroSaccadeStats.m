function setMicroSaccadeStats(obj)
    % Set the inter-saccde interval distribution
    obj.microSaccadeMeanIntervalSeconds = 0.45;
    obj.microSaccadeIntervalGammaShapeParameter = 5;
    
    % Set the microsaccade amplitude distribution 
    % (consistent with Fig 1 in Martinez 2009)
    obj.microSaccadeMeanAmplitudeArcMin = 8;
    obj.microSaccadeAmplitudeGammaShapeParameter = 4;
    
    % 30 deg/second (See Martinez 2008, Table in Fig 3D)
    obj.microSaccadeMeanSpeedDegsPerSecond = 39;
    obj.microSaccadeStDevSpeedDegsPerSecond = 2;  
    
    % Minimum duration of a microsaccade - Arbitrary value
    obj.microSaccadeMinDurationMilliSeconds = 2;
    
    % Sigma of microsaccade target jitter - Arbitrary value
    obj.microSaccadeTargetJitterArcMin = 0.3;
    
    % Sigma of corrective microsaccade direction jitter - Arbitrary value
    obj.microSaccadeDirectionJitterDegs = 15;
end