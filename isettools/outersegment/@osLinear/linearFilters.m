function [lmsFilters, meanCurrent] = linearFilters(os, cMosaic, varargin)
% Returns the photocurrent impulse response for a single absorption 
%
% The LMS impulse response functions calculated here model the cone
% photocurrent response to brief or low contrasts above and below a steady
% background.  These experimental conditions are often found in
% psychophysical or physiological experiment.
%
% The impulse response function is derived from Rieke's biophysical model.
% It depends on the mean absorption rate, and exhibits adaptation behavior.
% (s_osLinearFilters)
%
% There are different parameters for foveal and peripheral functions
% (t_osLinearize). 
% 
% Input
%  os:      A linear outer segment object
%  cMosaic: The parent object of the outersegment
%
% Return:
%  lmsFilters  - The impulse responses to a single photon added to
%                the background; these are stored in lmsConeFilter in
%                the os object
%  meanCurrent - The steady state caused by the mean absorption rate. This
%                value is used in osCompute(). 
%
% The LMS filters (impulse response functions) are stored here at a
% particular time step (os.timeStep), which is typically 1 ms, but could be
% shorter.  When it is used in osLinear.osCompute, the filters are
% resampled to the time base of the cone absorptions.
%
% See osLinear.osCompute() for how the impulse response function and mean
% currents are used.
%
% See also: v_osBioPhys, t_coneMosaicFoveal, t_osLinearize, s_osLinearFilters
%
% 11/2016 JRG/BW (c) isetbio team

%% parse input parameters
p = inputParser; p.KeepUnmatched = true;
p.addRequired('os', @(x) isa(x, 'outerSegment'));
p.addRequired('cMosaic', @(x) isa(x,'coneMosaic')); 

p.addParameter('eccentricity',false,@islogical);  % Needs updating - Foveal or peripheral
p.parse(os, cMosaic, varargin{:})

eccentricity   = p.Results.eccentricity;        % Needs updating
meanRate = coneMeanIsomerizations(cMosaic);     % R*/sec

%% Generate impulse responses for L, M or S cone

% A new cone mosaic is generated for each type of cone, and the
% impulse response
timeStep = os.timeStep;               % time step (should be < 1 ms)
nSamples = round(0.8/timeStep) + 1;   % 0.8 total sec

flashIntens = 1;   % 1 photon above the background mean
warmup = round(0.4/timeStep);    % Warm up period is 0.4 sec

% Where we store the filters
os.lmsConeFilter = zeros(nSamples-warmup+1,length(meanRate));
meanCurrent = zeros(1,3);

%% Generate a cone mosaic with an outerSegment based on the biophysical model 

% We turn off the noise and use the biophysical coneMosaic model to
% calculate an impulse response.
osCM = osBioPhys('osType',eccentricity);   % Will become eccentricity some day
osCM.set('noise flag','none');            % Run it without noise
cm = coneMosaic('os',osCM,'pattern', 2);   % single cone
cm.integrationTime = timeStep;
cm.os.timeStep = timeStep;

% For each of the cone types
for meanInd = 1:length(meanRate)
    
    % Get the isomerization rate (R*) in each time step R*
    meanIntens  = meanRate(meanInd)*timeStep;  
    
    % Create a constant stimulus at this rate
    stimulus = meanIntens*ones(nSamples, 1);
    
    % Compute outer segment current for the constant stimulus
    cm.absorptions  = reshape(stimulus, [1 1 nSamples]);
    cm.computeCurrent('bgR',meanRate(meanInd));
    currentConstant = squeeze(cm.current);
    % vcNewGraphWin; plot(timeStep*(1:nSamples),currentConstant);
    
    % Add a single photon (impulse) to the background one time step after
    % the warmup period
    stimulus(warmup+1) = meanIntens + flashIntens;
    
    % Compute outer segment currents with biophysical model
    % with the impulse stimulus
    cm.absorptions  = reshape(stimulus, [1 1 nSamples]);
    cm.computeCurrent('bgR',meanRate(meanInd));
    currentImpulse = squeeze(cm.current);
    % hold on; plot(timeStep*(1:nSamples),currentImpulse);
    
    % Store the impulse response.  We put flashIntens in for completeness,
    % but it is 1 so really, no need.
    os.lmsConeFilter(:,meanInd) = ...
        ((currentImpulse((warmup:end)-1)) - currentConstant((warmup:end)-1))/flashIntens;
    %vcNewGraphWin; 
    %plot(stimulus - meanIntens); hold on; 
    %plot(currentImpulse-currentConstant); grid on
    
    % We are a tiny bit worried about the edges; so we set a few steps
    % before the very end to the mean current from the constant stimulus
    % background. Edges are always a bitch.
    meanCurrent(meanInd) = currentConstant(end-10);
    
end

%% Assign them as output and return
lmsFilters = os.lmsConeFilter;

end