function em = emCreate(params)
%% Create eye movement structure
%
%   em = emCreate([params])
%
% Inputs:
%   params:  eye-movement parameters, optional, for default value see
%            auxilary function emFillParams below
%     .emType   - 3x1 flag vector, indicating whether or not to include
%                 tremor, drift and micro-saccade respectively
%     .sampTime - sampling time in secs, e.g. 0.001 stands for 1ms per
%                 sample
%     .totTime  - total time of eye-movement sequence in secs, will be
%                 rounded to a multiple of sampTime (removed)
%     .tremor   - parameters for tremor, could include
%        .interval   = mean of occurring interval in secs
%        .intervalSD = standard deviation of interval
%        .amplitude  = the randomness of tremor in rad
%        
%     .drift    - parameters for drift, could include
%        .speed     = speed of the slow drifting
%        .speedSD   = standard deviation of speed
%     .msaccade - parameters for microsaccade, could also be named as
%                 .microsaccade, could include fields as
%        .interval   = mean of occurring interval in secs
%        .intervalSD = standard deviation of interval
%        .dirSD      = the randomness of moving direction (rad)
%        .speed      = speed of micro-saccade (rad / sec)
%        .speedSD    = standard deviation of speed
%
% Output Parameter:
%   em       - eye movement structure
%
% Notes:
%   1) For all eye-movements, we assume that the acceleration time is
%      neglectable
%   2) Drift and tremor can be superimposed. However, during the 
%      microsaccade, both drift and tremor are suppressed
%   3) We assume that drift works like a 2D brownian motion. This is
%      reasonable when micro-saccade is present. However, when
%      micro-saccade is suppressed, drift will somehow towards the fixation
%      point. The actual pattern is unknown and this return-to-origin
%      feature is not implemented here.
%   4) For micro-saccade, speed and amplitude follows main-squence as
%      saccade and here, we just use speed
%   5) We don't add a field 'duration' to microsaccade and in computation,
%      the duration of microsaccade will be estimated by the speed and
%      distance between current position to the fixation point
%
% Reference:
%   1) Susana Martinez-Conde et. al, The role of fixational eye movements
%      in visual perception, Nature reviews | neuroscience, Vol. 5, 2004,
%      page 229~240
%   2) Susana Martinez-Conde et. al, Microsaccades: a neurophysiological
%      analysis, Trends in Neurosciences, Volume 32, Issue 9, September
%      2009, Pages 463~475
%
% Example:
%   em = emCreate;
%
% See also:
%   emGenSequence, emSet, emGet
%
% (HJ) Copyright PDCSOFT TEAM 2014

% Fill in default values for missing fields in params
if notDefined('params'), params = []; end

% set params to default values
% set general fields
p.name     = 'em structure';
p.type     = 'eye movement';
p.emFlag   = ones(3,1); % emType - no eye movement
p.sampTime = 0.001; % sample time interval - 1 ms

% set fields for tremor
p.tremor.interval   = 0.012;          % Tremor mean frequency - 83 Hz
p.tremor.intervalSD = 0.001;          % std of tremor frequency - 60~100 Hz
p.tremor.amplitude  = 18/3600*pi/180; % Tremor amplitude -  18 arcsec

% set fields for drift
% There's a big difference for drift speed between literatures, we just
% pick a reasonable value among them
p.drift.speed   = 3/60*pi/180;   % drift speed - drift speed, rad/sec
p.drift.speedSD = 1/60*pi/180;   % std of drift speed

% set fields for micro-saccades
p.msaccade.interval   = 0.6;       % micro-saccade interval - 0.6 secs
p.msaccade.intervalSD = 0.3;       % std for micro-saccade interval
p.msaccade.dirSD      = 5*pi/180;  % std for direction
p.msaccade.speed      = 15*pi/180; % micro saccade speed - 15 rad/s
p.msaccade.speedSD    = 5*pi/180;  % std for micro saccade speed

% merge params with default values
if ~isfield(params, 'msaccade') && isfield(params, 'microsaccade')
    params.msaccade = params.microsaccade;
    params = rmfield(params, 'microsaccade');
end
em = setstructfields(p, params);

% some checks for params
assert(numel(em.emFlag)==3, 'emType should be 3x1 logical vector');

end