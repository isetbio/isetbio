function varargout = v_osBioPhysEyeMovements(varargin)
% Validate the biophysical model of the cone outer segments against neural
% data obtained during a sequence of saccadic eye movements.
%
% This script tests the biophysically-based outer segment model of 
% photon isomerizations to photocurrent transduction that occurs in the
% cone outer segments.
%
% 6/xx/2015    npc   Created.
% 7/xx/2015    jrg   Test with ISETBIO outersegment object
% 12/31/15     dhb   Added local copy of coneAdapt.            
% 1/7/16       dhb   Rename.  Started to remove reference to coneAdapt.  
%                    Last version with coneAdapt comparison is in tagged
%                    version OSObjectVsOrigValidation.
% 1/12/16      npc   Created this version after separating the eye movements 
%                    component from s_coneModelValidate.

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

    %% Init
    ieInit;

    %% Load measured outer segment data
    [time, measuredOuterSegmentCurrent, stimulusPhotonRate] = loadMeasuredOuterSegmentResponses();
    
    %% Compute @os model response
    
    % Set the simulation time interval equal to the temporal sampling resolution of the measured measured data
    % In generar, the stimulation time interval should be set to a small enough value so as to avoid overflow errors.
    simulationTimeIntervalInSeconds = time(2)-time(1);
     
    % Create human sensor with 1 cone and load its photon rate with 
    % the stimulus photon rate time sequence
    sensor = sensorCreate('human');
    sensor = sensorSet(sensor, 'size', [1 1]); % only 1 cone
    sensor = sensorSet(sensor, 'time interval', simulationTimeIntervalInSeconds);
    sensor = sensorSet(sensor, 'photon rate', reshape(stimulusPhotonRate, [1 1 numel(stimulusPhotonRate)]));

    % Create a biophysically-based outersegment model object.
    osB = osBioPhys();
    
    % Specify no noise
    noiseFlag = 0;
    osB.osSet('noiseFlag', noiseFlag);

    % Compute the model's response to the stimulus
    osB.osCompute(sensor);

    % Get the computed current
    osBiophysOuterSegmentCurrent = osB.osGet('coneCurrentSignal');
    osBiophysOuterSegmentCurrent = squeeze(osBiophysOuterSegmentCurrent(1,1,:));
    
    offset1Time = 0.35;
    [~,offset1TimeBin] = min(abs(time - offset1Time ));

    offset2Time = 9.1;
    [~,offset2TimeBin] = min(abs(time - offset2Time ));
    
    % Add two different offset to measured current
    measuredOuterSegmentCurrentOffset1 = measuredOuterSegmentCurrent +  (osBiophysOuterSegmentCurrent(offset1TimeBin)-measuredOuterSegmentCurrent(offset1TimeBin));
    measuredOuterSegmentCurrentOffset2 = measuredOuterSegmentCurrent +  (osBiophysOuterSegmentCurrent(offset2TimeBin)-measuredOuterSegmentCurrent(offset2TimeBin));
    
    % compute RMS error
    residual1 = osBiophysOuterSegmentCurrent(:)-measuredOuterSegmentCurrentOffset1(:);
    residual2 = osBiophysOuterSegmentCurrent(:)-measuredOuterSegmentCurrentOffset2(:);
    validIndices = find(~isnan(measuredOuterSegmentCurrent));
    errorRMS1 = sqrt(mean(residual1(validIndices).^2));
    errorRMS2 = sqrt(mean(residual2(validIndices).^2));

    % Plot the two calculations and compare against measured data.
    if (runTimeParams.generatePlots)
        h = figure(1);
        set(h, 'Position', [10 1000 1000 1200]);
        subplot('Position', [0.05 0.54 0.94 0.42]);
        stairs(time,stimulusPhotonRate, 'r-',  'LineWidth', 2.0);
        set(gca, 'XLim', [time(1) time(end)], 'FontSize', 12);
        ylabel('Stimulus (R*/sec)','FontSize',14);
        
        subplot('Position', [0.05 0.03 0.94 0.46]);
        plot(time, measuredOuterSegmentCurrent, '.-', 'LineWidth', 2.0); hold on;
        plot(time, measuredOuterSegmentCurrentOffset1, 'm-', 'LineWidth', 2.0);
        plot(time, measuredOuterSegmentCurrentOffset2, 'b-', 'LineWidth', 2.0);
        plot(time, osBiophysOuterSegmentCurrent, 'k-',  'LineWidth', 2.0);
        plot(time(offset1TimeBin)*[1 1], [-100 100], 'm-');
        plot(time(offset2TimeBin)*[1 1], [-100 100], 'b-');
        set(gca, 'XLim', [time(1) time(end)], 'FontSize', 12);
        xlabel('Time (sec)','FontSize',14);
        ylabel('Photocurrent (pA)','FontSize',14);
        h = legend('measured (as saved in datafile)', sprintf('measured (adjusted to match model at %2.2f sec)', offset1Time),  sprintf('measured (adjusted to match model at %2.2f msec)',offset2Time) , 'osBioPhys model', 'location', 'NorthWest');
        set(h, 'FontSize', 12);
        title(sprintf('rms: %2.2f pA (offset at %2.2f sec)\nrms: %2.2f pA (offset at %2.2f sec)', errorRMS1, offset1Time, errorRMS2, offset2Time), 'FontName', 'Fixed');
        drawnow;
    end
    
    % Save validation data
    UnitTest.validationData('osBiophysCur', osBiophysOuterSegmentCurrent);
    UnitTest.validationData('time', time);
    UnitTest.validationData('stimulusPhotonRate', stimulusPhotonRate);
end

% Helper functions
function [time, measuredOuterSegmentCurrent, stimulusPhotonRate] = loadMeasuredOuterSegmentResponses()
    
    dataSource = {'resources/data/cones', 'eyeMovementExample'};
    fprintf('Fetching remote data: dir=''%s''  file=''%s''. Please wait ...\n', dataSource{1}, dataSource{2});
    % Download neural data from isetbio's repository
    client = RdtClient('isetbio');
    client.crp(dataSource{1});
    [eyeMovementExample, eyeMovementExampleArtifact] = client.readArtifact(dataSource{2}, 'type', 'mat');
    fprintf('Done fetching data.\n');
    
    extraTimeForBaselineComputation = 2.0;
    
    % time axis
    dt = eyeMovementExample.data.TimeAxis(2)-eyeMovementExample.data.TimeAxis(1);
    postStimulusTime = eyeMovementExample.data.TimeAxis(end) + dt*(1:(round(extraTimeForBaselineComputation/dt)));
    time = [eyeMovementExample.data.TimeAxis postStimulusTime];
    
    measuredOuterSegmentCurrent = nan(size(time));
    stimulusPhotonRate = time * 0;
    
    % Retrieve the (baseline-corrected) outer segment current
    stimTimeBins = 1:numel(eyeMovementExample.data.TimeAxis);
    measuredOuterSegmentCurrent(stimTimeBins) = squeeze(eyeMovementExample.data.Mean);
    
    % standard deviation of the current ?
    % measuredOuterSegmentCurrentSD = eyeMovementExample.data.SD;
    
    % stimulus in isomerizations/sec
    stimulusPhotonRate(stimTimeBins) = eyeMovementExample.data.Stim;
end

