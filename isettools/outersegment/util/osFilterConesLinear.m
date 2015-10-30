function newIRFs = osFilterConesLinear(varargin)
% Generates the linear temporal filters for the L, M and S cones based on
% physiological data from Angueyra and Rieke (Nat. Neuro., 2013). 
%
%   newIRFs = filterConesLinear([sensor]) 
%
% Inputs: sensor [optional] to match the sampling rate of the cone
% isomerization signal.
% 
% Outputs: 
%  newIRFs: nx3 matrix with each row the temporal impulse response for
%           the L, M and S cones, respectively.
%
% Examples:
% 
% newIRFs = osFilterConesLinear();
% newIRFs = osFilterConesLinear(sensor);
% 
% Originally by FMR
% Modified by JRG, isetbio team, 2015 


%% Handle input% 
% 
% TODO - Is it wrong to build the filters with a small number of samples?
% Does this matter for the kinds of simulation we are going to do?
% 
% Note from FMR - 
% Set the cone "sampling rate". This rate effectively dictates the rest of
% the simulation. Even though our monitor refreshed at 75 Hz, the cones
% were allowed to sample the stimulus more quickly. For convienence, I set
% the cone sampling rate to an integer multiple of the montor frame rate.
% Setting the cone sampling rate to 825 seemed reasonable given the spectra
% of the cone noise.

if isempty(varargin) || isempty(varargin{1})
    
    % Set default values when sensor is not given as input.
    dt = 0.001;
    coneSamplingRate = 1/dt;
    tsz = 300;
    meanIsomerization = 1;
    
else
    
    % Set values based on sensor properties
    sensor = varargin{1}; % sensor name
    dt = sensorGet(sensor, 'time interval');
    coneSamplingRate = 1/dt;
    [xsz ysz tsz] = size(sensor.data.volts);
    pRate = sensorGet(sensor, 'photons');
    meanIsomerization = mean(pRate(:));
    
end

%% Make the filter 
% Juan says the units are time (in sec) vs pA/R*/cone. All
% of the coefficients for this equation come from Juan's fits to
% electrophysiological measurements in a representative set of cones.

totalTime = dt * tsz; %length of IRF in seconds 

if totalTime <= 0.3; totalTime = 0.3; end;

TimeAxis= (0:ceil(coneSamplingRate.*totalTime)) ./ coneSamplingRate;


% only compute filter for 0.3 seconds
% TimeAxis = (1:timeBins)*dt;
TimeAxis = TimeAxis((TimeAxis <= 0.3));

% Parameters of cone temporal filters measured in Angueyra and Rieke (Nat.
% Neuro., 2013). These could conceivably be changed in the future, but not
% without further measurements from physiology.
ScFact = 1;      % 0.6745; % To get amplitude right
TauR = 0.0216;   % Rising Phase Time Constant
TauD = 0.0299;   % Damping Time Constant
TauP = 0.5311;   % Period
Phi = 34.1814;   % Phase
Filter = ScFact .* (((TimeAxis./TauR).^3)./(1+((TimeAxis./TauR).^3))) .* exp(-((TimeAxis./TauD))).*cos(((2.*pi.*TimeAxis)./TauP)+(2*pi*Phi/360));

%% Adjust gain of the cone impulse response function for L, M, S cones
%
% Note from FMR - 
% The IRF used by the model was calculated from real cones stimulated with
% white noise from a background of at ~13000 R*/sec. Adjust the gain of the
% model's IRF to reflect the adaptation state implied by the model's
% background light levels. By adjusting the gain, all I'm going to do is
% multiply the IRF by a scalar. This preserves the shape of the IRF across
% adaptation states, which isn't strictly empirically true but provides a
% reasonable approximation. There will be one gain factor for each cone
% type. Assume cone gain is a Weber-Fechner relationship with the
% half-desensitizing value indicated in Juan's paper (Angueyra & Rieke,
% Nat. Neuro., 2013). See Equation 1 of Juan's paper for details.
%
% Weber-Fechner relation:
%
%  gain / gain_dark = 1 / (1 + (Intensity_bkgnd / Intensitity_halfDesensitization))

Io = 2250;                     % half-desensitizing background (in R*/cone/sec, from Juan's paper - corrected)
Ib = [7131 6017 1973];         % R* per sec due to background adapting field (one for each cone, L, M, S)

stimNormCoeff = meanIsomerization/max(Ib);
% % %  If statement allows normalization to max of stimulus input:                              
% if size(varargin)==0
%     stimNormCoeff = meanIsomerization/max(Ib);
% else
% %     stimNormCoeff = (max(sensor.data.volts(:,:,1)))./max(Ib);%     
%     pRate = sensorGet(sensor,'photons');
%     stimNormCoeff = max(pRate(1))./max(Ib);
% end

% Set gain parameters according to Angueyra & Rieke:
Ib = Ib*stimNormCoeff;
gain_dark = 0.22;              % from Juan's paper (approximate peak of the IRF measured in darkness, and in units of pA/R*) - corrected
gainRatio = 1 ./ (1+(Ib./Io)); % the right side of the equation above, and the gain ratio implied by the bkgnd adapting field

% scale IRF to reflect amplitude at chosen background
% using Weber adaptation equation above and gainRatio derived from it
newGain = gainRatio .* gain_dark ;
oldGain = max(Filter);
IRFScaleFactor = newGain * dt ./ oldGain;

% The L, M and S IRFs are created by scaling Filter with the scale factors:
newIRFs = Filter(:) * IRFScaleFactor;

end
