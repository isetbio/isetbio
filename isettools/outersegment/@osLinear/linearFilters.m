function [lmsFilters, meanCurrent] = linearFilters(os, cMosaic, varargin)
% Returns the photocurrent impulse response for a single absorption 
%
% The impulse response function is calculated using the biophysical model.
% Hence, it depends on the mean absorption rate. 
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
% The impulse response functions calculated here are used when we are
% modeling a uniform background and stimuli that are smallish contrast
% modulations above and below that background.  These are often the
% conditions in a psychophysical or physiological experiment.
%
% See osLinear.osCompute() for how the impulse response function and mean
% currents are used.
%
% See also: v_osBioPhys, t_coneMosaicFoveal, t_osLinearize
%
% 11/2016 JRG/BW (c) isetbio team

% parse input parameters
p = inputParser; p.KeepUnmatched = true;
p.addRequired('os', @(x) isa(x, 'outerSegment'));
p.addRequired('cMosaic', @(x) isa(x,'coneMosaic')); 

p.addParameter('eccentricity',false,@islogical);  % Needs updating - Foveal or peripheral
p.parse(os, cMosaic, varargin{:})

eccentricity   = p.Results.eccentricity;        % Needs updating
meanRate = coneMeanIsomerizations(cMosaic);     % R*/sec

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

% Generate the cone mosaic with biophysical model for the impulse stimulus
osCM = osBioPhys('osType',eccentricity);   % Will become eccentricity some day
osCM.set('noise flag',0);                  % Run it without noise
cm = coneMosaic('os',osCM,'pattern', 2);   % single cone
cm.integrationTime = timeStep;
cm.os.timeStep = timeStep;

for meanInd = 1:length(meanRate)
    
    % Total isomerizations in each time step R* (maintained for 1 bin only)
    % We add one photon absorption to one mean
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
    
    % The mean current from the constant stimulus background.  We are a
    % tiny bit worried about the very end, so we take a few steps before
    % the very end.  Edges are always a bitch.
    meanCurrent(meanInd) = currentConstant(end-10);
    
end

lmsFilters = os.lmsConeFilter;

end