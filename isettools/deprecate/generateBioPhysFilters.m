function [lmsFilters, meanCurrent] = generateBioPhysFilters(os, meanRate, varargin)
% To be deprecated and replaced by linearFilters
%
% Returns the impulse response of the photocurrent given a single photon
% absorption 
%
% The impulse response function depends on the mean absorption rate.
%
% Input
%  os:         A linear outer segment object
%  meanRate:   a 3-vector with the photon absoprtion rates in the three
%              cone types (R*/sec)
%
% Return:
%  lmsFilters  - The impulse responses to a single photon added to
%                the background; these are also stored in lmsConeFilter in
%                the object itself
%  meanCurrent - The current in the steady state caused by the mean
%                 background
%
% The impulse response functions are used when we are modeling a uniform
% background and stimuli that are smallish contrast modulations above and
% below that background.  These are often the conditions in a
% psychophysical or physiological experiment.
%
% Manually sets an impulse stimulus in the absorptions field of a new cone
% mosaic object, generates the impulse response for the LMS cones and
% returns them, in order to be used in a linear computation with a
% near-constant background.
%
% See also: v_osBioPhys, t_coneMosaicFoveal, t_osLinearize
%
% 11/2016 JRG/BW (c) isetbio team

%%
disp('generateBioPhysFilters should be deprecated.  use linearFilters.')

%% parse input parameters
p = inputParser; p.KeepUnmatched = true;
p.addRequired('os', @(x) isa(x, 'outerSegment'));
p.addRequired('meanRate', @(x) length(x) == 3); % LMS mean absorption rates R*/sec
p.addParameter('osType',false,@islogical);  % Foveal or peripheral
p.parse(os, meanRate, varargin{:})

os       = p.Results.os;
osType   = p.Results.osType;
meanRate = p.Results.meanRate;

% Generate impulse responses for L, M or S cone
% A new cone mosaic is generated for each type of cone, and the
% impulse response
timeStep = os.timeStep;               % time step
nSamples = round(0.8/timeStep) + 1;   % 0.8 total sec

flashIntens = 1;   % 1 photon above the background mean
warmup = round(0.4/timeStep);    % Warm up period is 0.4 sec

% Where we store the filters
os.lmsConeFilter = zeros(nSamples-warmup+1,length(meanRate));
meanCurrent = zeros(1,3);

% Generate the cone mosaic with biophysical model for the
% impulse stimulus
osCM = osBioPhys('osType',osType);
osCM.set('noise flag','none');                  % Run it without noise
cm = coneMosaic('os',osCM,'pattern', 2);   % single cone
cm.integrationTime = timeStep;
cm.os.timeStep = timeStep;

for meanInd = 1:length(meanRate)
    
    % Total isomerizations in each time step R* (maintained for 1 bin only)
    meanIntens  = meanRate(meanInd)*timeStep;
    
    % Create a constant mean background stimulus.
    stimulus = meanIntens*ones(nSamples, 1);
    
    % Compute outer segment current for the constant stimulus
    cm.absorptions  = reshape(stimulus, [1 1 nSamples]);
    cm.computeCurrent('bgR',meanRate(meanInd));
    currentConstant = squeeze(cm.current);
    % vcNewGraphWin; plot(timeStep*(1:nSamples),currentConstant);
    
    % Add the impulse to the background one time step after the
    % warmup period
    stimulus(warmup+1) = meanIntens + flashIntens;
    
    % Compute outer segment currents with biophysical model
    % with the impulse stimulus
    cm.absorptions  = reshape(stimulus, [1 1 nSamples]);
    cm.computeCurrent('bgR',meanRate(meanInd));
    currentImpulse = squeeze(cm.current);
    % hold on; plot(timeStep*(1:nSamples),currentImpulse);
    
    % Check the +/-1 for warmup or warmup + 1 or ....
    os.lmsConeFilter(:,meanInd) = ...
        (currentImpulse((warmup:end)-1)) - currentConstant((warmup:end)-1) ./ flashIntens;
    % vcNewGraphWin; 
    % plot(stimulus - meanIntens); hold on; 
    % plot(currentImpulse-currentConstant); grid on
    % vcNewGraphWin; plot(impulseResponseLMS(:,meanInd)); grid on
    
    % The mean current from the constant stimulus background.  We are a
    % tiny bit worried about the very end, so we take a few steps before
    % the very end.  Edges are always a bitch.
    meanCurrent(meanInd) = currentConstant(end-10);
    
end

lmsFilters = os.lmsConeFilter;

end