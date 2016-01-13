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

    % remove DC offset
    osBiophysOuterSegmentCurrent = osBiophysOuterSegmentCurrent - osBiophysOuterSegmentCurrent(end);

    % compute RMS error
    residual = osBiophysOuterSegmentCurrent-measuredOuterSegmentCurrent;
    errorRMS = sqrt(mean(residual.^2));
    errorABS = mean(abs(residual));

    % Plot the two calculations and compare against measured data.
    if (runTimeParams.generatePlots)
        figure(1);
        subplot(2,1,1);
        stairs(time,stimulusPhotonRate, 'r-',  'LineWidth', 2.0);
        set(gca, 'XLim', [time(1) time(end)]);
        ylabel('Stimulus (photons/sec)','FontSize',14);
        
        subplot(2,1,2);
        plot(time, measuredOuterSegmentCurrent, 'm-', 'LineWidth', 2.0); hold on;
        plot(time, osBiophysOuterSegmentCurrent, 'k-',  'LineWidth', 2.0);
        set(gca, 'XLim', [time(1) time(end)]);
        xlabel('Time (sec)','FontSize',14);
        ylabel('Photocurrent (pA)','FontSize',14);
        legend('measured', 'osBioPhys model');
        title(sprintf('rms error: %2.2f pA', errorRMS));
        drawnow;
    end
    
    % Save validation data
    UnitTest.validationData('osBiophysCur', osBiophysOuterSegmentCurrent);
    UnitTest.validationData('time', time);
    UnitTest.validationData('stimulusPhotonRate', stimulusPhotonRate);
end

% Helper functions
function [time, measuredOuterSegmentCurrent, stimulusPhotonRate] = loadMeasuredOuterSegmentResponses()
    
    % Download neural data from isetbio's repository
    client = RdtClient('isetbio');
    client.crp('resources/data/cones');
    [eyeMovementExample, eyeMovementExampleArtifact] = client.readArtifact('eyeMovementExample', 'type', 'mat');

    % Retrieve the baseline corrected outer segment current
    measuredOuterSegmentCurrent = (squeeze(eyeMovementExample.data.Mean))';
    measuredOuterSegmentCurrent = measuredOuterSegmentCurrent - measuredOuterSegmentCurrent(end);
    
    % standard deviation of the current ?
    % measuredOuterSegmentCurrentSD = eyeMovementExample.data.SD;
    
    % stimulus in isomerizations/sec
    stimulusPhotonRate = eyeMovementExample.data.Stim;
    
    % time axis
    time = eyeMovementExample.data.TimeAxis;
end

