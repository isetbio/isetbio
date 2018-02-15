% Method to compute an emPath with duration emDurationSeconds
% sampled with a temporal resolution of sampleDurationSeconds
function computeSingleTrial(obj, emDurationSeconds, sampleDurationSeconds)
    
    % Init the state
    initState(obj, emDurationSeconds);
    
    feedbackDelaySteps(1) = obj.feedbackXposDelayTimeSteps;
    feedbackDelaySteps(2) = obj.feedbackYposDelayTimeSteps;

    % form time axis
    obj.timeAxis = ((1:obj.tStepsNum)-1) * obj.timeStepDurationSeconds;

    % Init the time step when the last saccade occurred
    obj.lastMicroSaccadeTimeStep = obj.stabilizationStepsNum;
         
    if (obj.displayComputeProgress)
        hWaitBar = waitbar(0, 'Computing emPath ...');
    end 
    
    for tStep = 1:obj.tStepsNum
        
        if (mod(tStep-1, 100) == 0) && (obj.displayComputeProgress)
            progress = tStep/obj.tStepsNum;
            waitbar(progress, hWaitBar, 'Computing emPath ...');
        end
        for eyeIndex = 1:2
            % compute feedback signals
            if (tStep > feedbackDelaySteps)
                delayedControl = obj.controlSignalTimeSeries(eyeIndex, tStep-feedbackDelaySteps(eyeIndex));
                obj.feedbackSignalTimeSeries(eyeIndex,tStep) = ...
                    obj.feedbackGain * tanh(obj.feedbackSteepness*delayedControl);
            end

            % compute next control signals
            obj.controlSignalTimeSeries(eyeIndex,tStep+1) = ...
                (1-obj.controlGamma)*obj.controlSignalTimeSeries(eyeIndex,tStep) - ...
                obj.feedbackSignalTimeSeries(eyeIndex,tStep) + obj.controlNoiseTimeSeries(eyeIndex,tStep);

            % generate next EMpoint
            obj.emPosTimeSeries(eyeIndex,tStep+1) = ...
                obj.emPosTimeSeries(eyeIndex,tStep) + ...
                obj.controlSignalTimeSeries(eyeIndex,tStep+1) + ...
                obj.positionalNoiseTimeSeries(eyeIndex, tStep) + ...
                obj.microSaccadeResidualPath(eyeIndex,1);
        end % eyeIndex

        % Update the microSaccadeResidualPath
        obj.updateMicroSaccadeResidualPath();
            
        % Check for saccade generation only after the drift has stabilized
        if (tStep>=obj.stabilizationStepsNum)
            if strcmp(obj.microSaccadeType ,'heatmap/fixation based')
                % Update heat map and retrieve heat level at the current emPos. 
                % Maybe we could also trigger saccades when the
                % currentPositionHeatLevel exceeds a threshold.
                currentPositionHeatLevel = obj.updateHeatMap(tStep);
            end
            if (~strcmp(obj.microSaccadeType,'none'))
                % Check whether a microsaccade will be generated
                obj.checkForMicroSaccadeEpoch(tStep);
            end
        end % if (tStep>=obj.stabilizationStepsNum)
    end % tStep

    % Compute velocity
    obj.velocityTimeSeries = obj.computeVelocity(obj.emPosTimeSeries);
    
    % Trim and recenter and resample the post-stabilization time series
    trimRecenterAndResampleTimeSeries(obj, sampleDurationSeconds);

    if (strcmp(obj.microSaccadeType, 'heatmap/fixation based'))
        % Scale heat map to 1 (this affects only visualization at this point)
        obj.heatMapTimeSeries = obj.heatMapTimeSeries / max(obj.heatMapTimeSeries(:));
    end
    
    % Make the arc min sequence
    obj.emPosTimeSeriesArcMin = obj.emPosTimeSeries*obj.scalarToArcMin;   
    obj.velocityArcMinPerSecTimeSeries = obj.velocityTimeSeries*obj.scalarToArcMin;
    
    if (obj.displayComputeProgress) && (obj.displayComputeProgress)
        close(hWaitBar);
    end
end % compute 