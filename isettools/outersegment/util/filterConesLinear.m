function [newIRFs, Filter, Ib] = filterConesLinear(varargin)


%%  CODE TO GENERATE THE CONE TEMPORAL IMPULSE RESPONSE FUNCTION
%
% Set the cone "sampling rate". This rate effectively dictates the rest of
% the simulation. Even though our monitor refreshed at 75 Hz, the cones
% were allowed to sample the stimulus more quickly. For convienence, I set
% the cone sampling rate to an integer multiple of the montor frame rate.
% Setting the cone sampling rate to 825 seemed reasonable given the spectra
% of the cone noise.

if size(varargin)==0
    dt = 0.001;
else
    sensor = varargin{1}; % sensor name
    dt = sensorGet(sensor, 'time interval');
end

%make the filter. Juan says the units are time (in sec) vs pA/R*/cone. All
%of the coefficients for this equation come from Juan's fits to
%electrophysiological measurements in a representative set of cones.
totalTime = 2; %length of IRF in seconds 
% TimeAxis= (0:ceil(coneSamplingRate.*totalTime)) ./ coneSamplingRate;
TimeAxis = [1:2000]*dt;
ScFact = 0.6745; % To get amplitude right
TauR = 0.0216;   % Rising Phase Time Constant
TauD = 0.0299;   % Damping Time Constant
TauP = 0.5311;   % Period
Phi = 34.1814;   % Phase
Filter = ScFact .* (((TimeAxis./TauR).^3)./(1+((TimeAxis./TauR).^3))) .* exp(-((TimeAxis./TauD))).*cos(((2.*pi.*TimeAxis)./TauP)+(2*pi*Phi/360));

%% CODE TO ADJUST GAIN OF THE CONE IRF
%
% The IRF used by the model was calculated from real cones stimulated with
% white noise from a background of at ~13000 R*/sec. Adjust the gain of the
% model's IRF to reflect the adaptation state implied by the model's
% background light levels. By adjusting the gain, all I'm going to do is
% multiply the IRF by a scalar. This preserves the shape of the IRF across
% adaptation states, which isn't strictly empirically true but provides a
% reasonable approximation. There will be one gain factor for each cone
% type. Assume cone gain is a Weber-Fechner relationship with the
% half-desensitizing value indicated in Juan's paper (Angueyra & Rieke
% 2013). See Equation 1 of Juan's paper for details.
%
% Weber-Fechner relation:
%
%  gain / gain_dark = 1 / (1 + (Intensity_bkgnd / Intensitity_halfDesensitization))

Io = 2250;                     % half-desensitizing background (in R*/cone/sec, from Juan's paper - corrected)
Ib = [7131 6017 1973];         % R* per sec due to background adapting field (one for each cone, L, M, S)
% Ib = [7131 7131 7131];  
                                % adjust this to specific experiment
if size(varargin)==0
    stimNormCoeff = 1;
else
%     stimNormCoeff = (max(sensor.data.volts(:,:,1)))./max(Ib);
    
    pRate = sensorGet(sensor,'photon rate');
    stimNormCoeff = max(pRate(:))./max(Ib);
end

Ib = Ib*stimNormCoeff;
gain_dark = 0.32;              % from Juan's paper (approximate peak of the IRF measured in darkness, and in units of pA/R*) - corrected
gainRatio = 1 ./ (1+(Ib./Io)); % the right side of the equation above, and the gain ratio implied by the bkgnd adapting field



% scale IRF to reflect amplitude at chosen background
% using Weber adaptation equation above and gainRatio derived from it
newGain = gainRatio .* gain_dark ;
oldGain = max(Filter);
IRFScaleFactor = newGain * dt ./ oldGain;


% plot the original IRF, and the IRFs that have been adjusted to reflect
% their particular adaptation state
newIRFs = Filter(:) * IRFScaleFactor;

% 
% figure(1); clf
% plot(TimeAxis, newIRFs(:, [3,2,1]))  % implied IRFs at the background intensities used in the model
% set(gca, 'fontsize', 14)
% legend('IRF S-Cones', 'IRF M-Cones', 'IRF L-Cones')
% xlabel('Time (sec)')
% ylabel('pA per R*')
% title('Cone IRFs at specified background')
% xlim([0 0.2]);

%%
% convolve with sinusoidal stimulus to generate signal

Frequency = 20;          % Hz
Contrast = 0.1;         
Ib = 7131;              % L cone
% Filt = newIRFs(:, 1)';
% TaperLength = round(length(Filt) / 4);
% 
% StimEnvelope = 1:length(Filt);
% StimEnvelope(:) = 1;
% StimEnvelope(1:TaperLength) = (0:TaperLength-1)*1/TaperLength;
% StimEnvelope(length(Filt)-TaperLength+1:length(Filt)) = 1-(0:TaperLength-1)*1/TaperLength;
% 
% Stim = Contrast * Ib * sin(2*3.1416*Frequency*TimeAxis);
% Stim = Stim .* StimEnvelope;
% 
% % convolve
% ConeSignal = real(ifft(fft(Filt) .* fft(Stim))) / coneSamplingRate; 
% 
% figure(3);
% plot(TimeAxis, ConeSignal);
% xlabel('time (sec)');
% ylabel('pA');
% 
% %%
% % add noise
% 
% FreqStep = coneSamplingRate / length(ConeSignal);
% Freq = (1:length(ConeSignal)) * FreqStep;
% 
% LorentzCoeffs=[0.16   55    4    0.045  190    2.5];
% Fit=lorentzsum_poles(LorentzCoeffs,Freq);
% 
% % generate white gaussian noise of correct length
% Noise = normrnd(0, 1, size(ConeSignal));
% ConeNoiseFFT = fft(Noise);
% 
% for i=1:length(Noise)/2
%     ConeNoiseFFT(i) = ConeNoiseFFT(i) * sqrt(Fit(i));
%     ConeNoiseFFT(length(Noise) - i + 1) = ConeNoiseFFT(length(Noise) - i + 1) * sqrt(Fit(i));
% end
% 
% ConeNoise = real(ifft(ConeNoiseFFT)) * sqrt(coneSamplingRate) / sqrt(2); 
% 
% figure(4);
% loglog(Freq, Fit, Freq, 2 * real(fft(ConeNoise) .* conj(fft(ConeNoise))) / (coneSamplingRate * length(ConeNoise)));