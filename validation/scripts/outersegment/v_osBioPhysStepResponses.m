function varargout = v_osBioPhysStepResponses(varargin)
% Validate the biophysical model of the cone outer segments against neural
% data obtained in response to step stimuli of different amplitudes.
%
% This script tests the biophysically-based outer segment model of 
% photon isomerizations to photocurrent transduction that occurs in the
% cone outer segments.
%
% 1/12/16      npc   Created after separating the relevant 
%                    components from s_coneAdaptNoise.

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

    %% Init
    ieInit;

    %% Load measured outer segment data
    [time, measuredOuterSegmentCurrents, stimulusPhotonRates] = loadMeasuredOuterSegmentResponses();
    
    % Set the simulation time interval equal to the temporal sampling resolution of the measured measured data
    % In general, the stimulation time interval should be set to a small enough value so as to avoid overflow errors.
    upsamplingFactor = 1;
    simulationTimeIntervalInSeconds = (time(2)-time(1))/upsamplingFactor;
    simulationTime = time(1):simulationTimeIntervalInSeconds:time(end);
    nSamples = numel(simulationTime);
    pulseOnset  =  40001;
    pulseOffset = 120000;
    
    % create human sensor with 1 cone
    sensor = sensorCreate('human');
    sensor = sensorSet(sensor, 'size', [1 1]); % only 1 cone
    sensor = sensorSet(sensor, 'time interval', simulationTimeIntervalInSeconds);
        
    for stepIndex = 1:numel(stimulusPhotonRates)
        
        % retrieve the measured currents for the examined step size
        measuredCurrents = measuredOuterSegmentCurrents{stepIndex};
        
        % retrieve the examined stimulus amplitude (photons/sec)
        stimulusPhotonRateAmplitude = stimulusPhotonRates{stepIndex};
        
        % start and end time of the light
        stimPeriod = [pulseOnset pulseOffset]*upsamplingFactor; 

        % create stimulus temporal profile
        stimulusPhotonRate = zeros(nSamples, 1);
        stimulusPhotonRate(stimPeriod(1):stimPeriod(2)) = stimulusPhotonRateAmplitude;
        
        % set the stimulus photon rate
        sensor = sensorSet(sensor, 'photon rate', reshape(stimulusPhotonRate, [1 1 numel(stimulusPhotonRate)]));

        % create a biophysically-based outersegment model object
        osB = osBioPhys();
        
        % specify no noise
        noiseFlag = 0;
        osB.osSet('noiseFlag', noiseFlag);
    
        % compute the model's response to the stimulus
        osB.osCompute(sensor);
    
        % get the computed current
        current = osB.osGet('coneCurrentSignal');
        
        % store copy for saving to validation file
        if (stepIndex == 1)
            osBiophysOuterSegmentCurrent = zeros(numel(stimulusPhotonRates), size(current,3));
        end
        osBiophysOuterSegmentCurrent(stepIndex,:) = squeeze(current(1,1,:));
    
        % plot model and measured data and stimulus
        if (runTimeParams.generatePlots)
            if (stepIndex == 1)
                h = figure(1); clf;
                set(h, 'Position', [10 10 900 1200]);
            end
            
            % plot stimulus on the left
            subplot(numel(stimulusPhotonRates),2,(stepIndex-1)*2+1);
            stairs(simulationTime, stimulusPhotonRate, 'r-', 'LineWidth', 2.0);
            set(gca, 'YLim', [0 12e5]);
            if (stepIndex == numel(stimulusPhotonRates))
            	xlabel('time (sec)','FontSize',12);
            else
                set(gca, 'XTickLabel', {});
            end
           
            ylabel('stimulus (R*/sec)','FontSize',12);
            text(2.4, 9e5, sprintf('step=%d R*/sec',stimulusPhotonRateAmplitude), 'FontSize',12);
        
            % plot responses on the right
            subplot(numel(stimulusPhotonRates),2,(stepIndex-1)*2+2);
            hold on;
            trialColors = [...
                1.0 0.6 0.2; ...
                0.6 1.0 0.5; ...
                0.3 0.4 1.0; ...
                0.9 0.2 1.0 ...
                ];
            % plot individual trial measured  currents offset
            for trial = 1:size(measuredCurrents,1)
                measuredOuterSegmentCurrent = squeeze(measuredCurrents(trial,:));
                plot(time, measuredOuterSegmentCurrent, '-', 'Color', squeeze(trialColors(trial,:)), 'LineWidth', 2.0);
            end
            
            % plot computed current
            plot(simulationTime, squeeze(osBiophysOuterSegmentCurrent(stepIndex,:)), 'k-', 'LineWidth', 2.0);
            if (stepIndex == numel(stimulusPhotonRates))
            	xlabel('time (sec)','FontSize',12);
            else
                set(gca, 'XTickLabel', {});
            end
            
            %set(gca, 'YLim', [-150 0]);
            ylabel('current (pAmps)','FontSize',12);
            if (size(measuredCurrents,1) == 1)
                legend('measured (trial-1)', 'osBioPhys model');
            elseif (size(measuredCurrents,1) == 2)
                legend('measured (trial-1)', 'measured (trial-2)', 'osBioPhys model');
            elseif (size(measuredCurrents,1) == 3)
                legend('measured (trial-1)', 'measured (trial-2)', 'measured (trial-3)', 'osBioPhys model');
            end
            
%             if (stepIndex == numel(stimulusPhotonRates))
%                 NicePlot.exportFigToPNG('PulseTests2.png', h, 300);
%             end
        end   
    end % stepIndex
    
    % Save validation data
    UnitTest.validationData('osBiophysCurrent', osBiophysOuterSegmentCurrent);
    UnitTest.validationData('simulationTime', simulationTime);
    UnitTest.validationData('stimPeriod', stimPeriod);
    UnitTest.validationData('stimulusPhotonRates', stimulusPhotonRates);
end

% Helper functions
function [time, measuredOuterSegmentCurrents, stimulusPhotonRates] = loadMeasuredOuterSegmentResponses()
    
    dataSource = {'resources/data/cones', 'stepExample'};
    fprintf('Fetching remote data: dir=''%s''  file=''%s''. Please wait ...\n', dataSource{1}, dataSource{2});
    % Download neural data from isetbio's repository
    client = RdtClient('isetbio');
    client.crp(dataSource{1});
    [stepExample, stepExampleArtifact] = client.readArtifact(dataSource{2}, 'type', 'mat');
    fprintf('Done fetching data.\n');

    % stimulus in isomerizations/sec
    stimulusPhotonRates = stepExample.data.lightLevel;
    
    
    % measured outer segment currents
    measuredOuterSegmentCurrents = stepExample.data.Data;
    
    % time axis
    time = (1:size(measuredOuterSegmentCurrents{1},2))*stepExample.data.samplingInterval;
end

