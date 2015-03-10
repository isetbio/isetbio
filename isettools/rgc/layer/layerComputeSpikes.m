function [spkTimes, rgcV, rgcThresh] = layerComputeSpikes(layer)
% Compute spikes for a layer
% 
%  [spkTimes, rgcV, rgcThresh] = layerComputeSpikes(layer)
%
% Returns:
%  - spkTimes:  The spike time series
%  - rgcV:      Matrix of RGC voltages that include input, noise, feedback.
%  - rgcThresh: Voltage level for spiking event.
%
%  These are updated at every time step.
%
% Algorithm:
%  0. The RGC voltage (rgcV) is the sum of the linear time series
%  (linTS), the feedback signals caused by spikes (spkTS), and the noise
%  time series (noiseTS).
%
%  1. If the combined signal exceeds a threshold, a spike is created.  
%
%  2. Spikes initiate two feedback events:
%    A. fbTR (auto feedback time response, vector): 
%      This voltage, which extends over time, is added to that RGC's
%      membrane potential. Each time a neuron spikes, we start a new
%      feedback response voltage and add it to spkTS for that RGC, starting
%      at the current time step. 
%    B. cpTR (coupling signal, vector):
%     This coupling voltage is sent to connected neurons. The connected
%     RGCS are specified by the connection matrix. The cpTR is added to
%     neighboring neurons' spkTS.
%
% At each time step, linTS, spkTS and noiseTS are added to the current
% membrane potential rgcV to make the rgcV of the next time step. 
%
% The connection matrix (connections) contains, for each neuron, the list
% of its connected neighbors, and the gain of that connection. The farther
% away the neighbor is, the lower the gain is. 
%
% The feedback impulse response (fbTR) issued by a spiking neuron is
% multiplied by this gain before being added to the neighbor neuron's
% spkTS. 
%
% See wiki page (may be obsolete by now)
%  http://white.stanford.edu/pdcwiki/index.php/Spike_generation_using_a_Poisson_process
%
% Example:
%
% See also: rgcComputeSpikes
%
% (c) 2010 Stanford Vistasoft Team

% Programming TODO:
%   The spikes are saturating.

%% Input checking
if notDefined('layer'), error('RGC layer structure required'); end

% We will need parameters from the parent RGC structure
rgcP        = layer.parent;
connections = layer.get('Connection Matrix');

% The linear time series has a minimum of 0 and runs into the general size
% of the RGC voltage swing.
linTS       = layer.get('Linear Time Series');

% Should we rectify the linear time series?  If so, should there be a
% baseline?  This should be handled in rgcAdapt.

% Set the voltage level (rgcThresh) for a spike
% RGC Volt Thresh is a fraction that we apply to the voltage swing.
rgcThresh = layer.get('RGC Volt Thresh');  

%% Retrieve layer and rgcP parameters

% Temporal parameters
stimTime    = rgcP.get('nFrames');         % Number of timepoints

% Retrieve network parameters
gridSize = layer.get('gridSize');  % Sample size of the RGCs in this grid - Check this more.
nRGC     = prod(gridSize);         % Number of RGCs

% Coupling and feedback. If these properties are turned off by
%  layer.disableCoupling or layer.hasFeedback = 0, then these should be 0
%  and thus have no effect.
cpTR     = layer.get('cpTR');      % Coupling signal temporal response
fbTR     = layer.get('fbTR');      % Feedback impulse response
lastFeedbackTime = max(length(cpTR),length(fbTR));   % Longest feedback

% Debugging
% figure(1); plot(fbTR); grid on; xlabel('Time (ms)'); ylabel('Resp (mv)')
% figure(1); plot(cpTR); grid on
%
% Visualize linTS if you like ...
% Not sure how to generalize.  We need row and col.  Then we get a video of
% the linear time series kind of like this.  Or through the other routine
% that we have for visualizing videos.
% for kk=1:size(linTS,2)
%     tmp = linTS(:,kk);
%     imagesc(reshape(tmp,gridSize(1),gridSize(2)))
%     pause(0.1)
% end

%% Initialize rgcV time series, noiseTS, spike times, coupling

% Retrieve noise parameters and create noise time series
% mV  = rgcP.get('meanV');  % Mean voltage at the RGC or Cone?
% sV  = rgcP.get('stdV');   % Std voltage 

% This is Gaussian noise to be added to the RGC voltage
% Commented out below.  Perhaps because of data from Rieke.
% dT    = rgcP.get('dT');              % Simulation timestep ms
% noiseTS  = ( (randn(size(linTS)) * sV) + mV) * dT;  
   
% Remembered integral state from preceding simulation (size rGC*1)
% spikeIntegral = layer.get('currentSpikeIntegral');
% Probably need to delete this from the layer structure.

% We store the voltages caused by feedback temporal response (fbTR) and network
% coupling temporal response (cpTR) in these time series.
% These time series are initialized to 0
% fbTRTS = zeros(nRGC, stimTime+lastFeedbackTime);
% cpTRTS = zeros(nRGC, stimTime+lastFeedbackTime);

totalTime = stimTime + lastFeedbackTime - 1;

% Matrix representing spike times.  0s and 1s.
spkTimes = zeros(nRGC,totalTime,'uint8');

% Matrix representing voltage created from spikes (feedback and coupling)
spkTS = zeros(nRGC,totalTime);

% Initialize the rgc time series as the linear time series.
rgcV = linTS;

%% Run simulation time-step by time-step

% We count from t+1 below, so we can only go up to stimTIme - 1.
% This is a slight temporal edge we should worry about.
for t = 1:(stimTime-1)   

    % There are alternative models for how the RGC gets updated.  We want
    % the code to reflect some sluggishness in the cones, and separately a
    % temporal integration on the RGC integrated signal as well.  So we
    % should experiment with models like this ...
    % rgcV(:,t) = rgcV(:,t-1)*leak + linTS(:,t) + spkTS(:, t) + noiseTS(:,t);
    
    % NOTE:  
    % We want to add the noise to the linTS so we don't disturb the SNR of
    % the process.  After adding the noise to the signal, we should scale
    % (contrast adaptation) so that the total signal is in the same
    % operating range that we expect to work with the spike-initiated
    % feedback properties.
    
    % But for now we are using the following, which is not good.
    % For debugging, we eliminate the noise.  There is some chance we will
    % always eliminate it because of the Rieke-Chichilnisky paper.
    % rgcV(:,t) = linTS(:,t) + spkTS(:, t) + noiseTS(:,t);     
    rgcV(:,t) = rgcV(:,t) + spkTS(:, t);     
    % figure; plot(linTS(:,t))
    
    % Across the RGCs, determine the voltage level: Is this a spike?
    % We have to deal with positive and negative voltage somewhere.
    spkTimes(:,t) = (rgcV(:,t) > rgcThresh); 
    % max(rgcV(:,t));  figure; hist(rgcV(:,t),50)
    
    % Find RGCs that have spiked 
    idxRGC  = (spkTimes(:,t) == 1);
    
    % lst is the numbers of the RGCs where there was a spike.
    lst = find(idxRGC(:) ~= 0);
    
    % If no spikes, break out of loop
    if isequal(idxRGC,zeros(nRGC,1)), continue; end;  % 
    
    % At least one spike occurred, so add the feedback and coupling signals
    % to the spkTS variable (if they are enabled).  This is the spike
    % driven time series.  
    
    if layer.hasFeedback
        % Calculate the indices of future timepoints where we add feedback.
        iTR   = (t+1) : (t + length(fbTR));

        % Add this layer's feedback function, fbTR, to the spkTS at the future
        % timepoints
        spkTS(idxRGC, iTR)  = spkTS(idxRGC, iTR)  + repmat(fbTR, sum(idxRGC),1);
        % figure; plot(fbTR); figure; imagesc(spkTS)

    end
    
    
    if layer.hasCoupling
        % Add the coupling signals to spkTS for those neurons connected to the
        % spiking neurons.

        % Because the length of the coupling time course may differ, we repeat
        % the whole calculation.  If we force them to have the same length, we
        % could skip (only) this line.
        iTR      = (t + 1):(t + length(cpTR));
        
        % Build the cpSignal (coupling) matrix to add to spkTS
        % Initialize the coupling signal.  And then loop through the RGCs that
        % spiked.
        cpSignal = zeros(nRGC,length(iTR));
        % figure; plot(cpTR)
        
        for ii = 1:length(lst)
            % For every RGC, find its connections and add the coupling signal
            % to those connected RGCs.
            if ~isempty(connections{lst(ii),1})  % If it is connected
                jj = connections{lst(ii),1};     % Find its connections
                % The second entry of connections is the weighted cpTR between
                % the lst(ii) neuron and the jj neuron.
                cpSignal(jj,:) = cpSignal(jj,:) + connections{lst(ii),2};
            end
        end

        % Add coupling signal into the spike time series
        spkTS(:, iTR) = spkTS(:, iTR) + cpSignal;
    end

end

% updating current state:
% The combined linTS, spkTS, and noiseTS
layer.set('rgc voltage time series',rgcV);

% The spike-initiated time series
layer.set('currentSpkTS',spkTS);

return
