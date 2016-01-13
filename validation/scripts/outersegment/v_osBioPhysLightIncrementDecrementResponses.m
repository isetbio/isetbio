function varargout = v_osBioPhysLightIncrementDecrementResponses(varargin)
% Validate the biophysical model for light increment and decrement stimuli
%
% This script tests the biophysically-based outer segment model of 
% photon isomerizations to photocurrent transduction that occurs in the
% cone outer segments.
%
% 1/12/16      npc   Created after separating the relevant 
%                    components from s_coneModelValidate.

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

    %% Init
    ieInit;
    
    % Set the simulation time interval. In general, the stimulation time interval should 
    % be set to a small enough value so as to avoid overflow errors.
    simulationTimeIntervalInSeconds = 5e-4;
    
    % Compute the simulation time axis
    pulseOnset  = 4000;
    pulseOffset = 7500;
    stimPeriod = [pulseOnset pulseOffset];
    nSamples   = pulseOffset+2000;
    simulationTime = (1:nSamples)*simulationTimeIntervalInSeconds;
    
    stimulusPhotonRateAmplitudes = 500 * 2.^(1:7); % photons/sec
    contrastsExamined = [-1 1];
    
    for stepIndex = 1:numel(stimulusPhotonRateAmplitudes)
        
        % create stimulus temporal profile
        stimulusPhotonRate = zeros(nSamples, 1);
        stimulusPhotonRate(100:nSamples-100,1) = stimulusPhotonRateAmplitudes(stepIndex);
        
        for contrastIndex = 1:numel(contrastsExamined)   

            fprintf('Running simulation for step #%d, contrast #: %d\n', stepIndex, contrastIndex);
            
            % generate step (decrement/increment)
            stimulusPhotonRateStep(contrastIndex, :) = stimulusPhotonRate;
            stimulusPhotonRateStep(contrastIndex, pulseOnset:pulseOffset) = stimulusPhotonRate(pulseOnset:pulseOffset,1) * (1+contrastsExamined(contrastIndex));
            
            % create human sensor with 1 cone and load its photon rate with 
            % the stimulus photon rate time sequence
            sensor = sensorCreate('human');
            sensor = sensorSet(sensor, 'size', [1 1]); % only 1 cone
            sensor = sensorSet(sensor, 'time interval', simulationTimeIntervalInSeconds);
            
            % set the stimulus photon rate
            sensor = sensorSet(sensor, 'photon rate', reshape(squeeze(stimulusPhotonRateStep(contrastIndex,:)), [1 1 size(stimulusPhotonRateStep,2)]));

            
            % create a biophysically-based outersegment model object
            osB = osBioPhys();
        
            % specify no noise
            noiseFlag = 0;
            osB.osSet('noiseFlag', noiseFlag);
    
            % compute the model's response to the stimulus
            params.bgVolts = 0;
            osB.osCompute(sensor, params);
            
            % get the computed current
            current = osB.osGet('coneCurrentSignal');

            % store copy for saving to validation file
            if ((stepIndex == 1) && (contrastIndex == 1))
                osBiophysOuterSegmentCurrent = zeros(numel(stimulusPhotonRateAmplitudes), numel(contrastsExamined), size(current,3));
            end
            osBiophysOuterSegmentCurrent(stepIndex, contrastIndex,:) = current(1,1,:);
            
        end % contrastIndex
        
        if (runTimeParams.generatePlots)
            
            if ((stepIndex == 1))
                h = figure(1); clf;
                set(h, 'Position', [10 10 900 1200]);
            end
            
            % plot stimulus on the left
            indices = 1:2:numel(simulationTime);
            subplot(numel(stimulusPhotonRateAmplitudes),2,(stepIndex-1)*2+1); hold on;
            plot([simulationTime(1) simulationTime(end)], stimulusPhotonRate*[1 1], 'k-');
            stairs(simulationTime(indices), squeeze(stimulusPhotonRateStep(1,indices)), 'm-', 'LineWidth', 2.0);
            stairs(simulationTime(indices), squeeze(stimulusPhotonRateStep(2,indices)), 'b-', 'LineWidth', 1.0);
            set(gca, 'XLim', [simulationTime(1) simulationTime(end)], 'YLim', [0 15e4]);
            if (stepIndex == numel(stimulusPhotonRateAmplitudes))
                xlabel('time (sec)','FontSize',12);
            else
                set(gca, 'XTickLabel', {});
            end
            ylabel('stimulus (photons/sec)','FontSize',12);
            text(0.1, 10e4, sprintf('pedestal: %d photons/sec',stimulusPhotonRateAmplitudes(stepIndex)), 'FontSize',12);

            % plot responses on the right
            subplot(numel(stimulusPhotonRateAmplitudes),2,(stepIndex-1)*2+2); hold on;
            plot(simulationTime(indices), squeeze(osBiophysOuterSegmentCurrent(stepIndex, 1, indices)), 'm-', 'LineWidth', 2.0);
            plot(simulationTime(indices), squeeze(osBiophysOuterSegmentCurrent(stepIndex, 2, indices)), 'b-', 'LineWidth', 1.0);
            set(gca, 'XLim', [simulationTime(1) simulationTime(end)], 'YLim', [-100 0]);
            if (stepIndex == numel(stimulusPhotonRateAmplitudes))
                xlabel('time (sec)','FontSize',12);
            else
                set(gca, 'XTickLabel', {});
            end
            ylabel('current (uAmps)','FontSize',12);
        end % if (runTimeParams.generatePlots)
    end % stepIndex
    
    fprintf('Generating figure. This will take a while ...\n');
    tic
    drawnow;
    fprintf('Figure generation took %2.2f sec\n', toc);
    
    % Save validation data
    UnitTest.validationData('osBiophysCurrent', osBiophysOuterSegmentCurrent);
    UnitTest.validationData('simulationTime', simulationTime);
    UnitTest.validationData('stimPeriod', stimPeriod);
    UnitTest.validationData('stimulusPhotonRateAmplitudes',stimulusPhotonRateAmplitudes);
end

