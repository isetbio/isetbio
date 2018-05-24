function em = emCreate(params, varargin)
% Create eye movement structure
%
% Syntax:
%   em = emCreate(params, [varargin])
%
% Description:
%    Create the structure to contain the eye movement.
%
% Inputs:
%    params   - eye-movement parameters
%      .emFlag   - 3x1 flag vector, indicating whether or not to include
%                  tremor, drift and micro-saccade respectively
%      .sampTime - sampling time in secs, e.g. 0.001 stands for 1ms per
%                  sample
%      .tremor   - parameters for tremor, could include
%         .interval = mean of occurring interval in secs
%         .intervalSD = standard deviation of interval
%         .amplitude = the randomness of tremor in radians
%      .drift    - parameters for drift, could include
%         .speed = speed of the slow drifting
%         .speedSD = standard deviation of speed
%      .msaccade - parameters for microsaccade, could also be named as
%                  .microsaccade, could include fields as
%         .interval = mean of occurring interval in secs
%         .intervalSD = standard deviation of interval
%         .dirSD = the randomness of moving direction (rad)
%         .speed = speed of micro-saccade (rad / sec)
%         .speedSD = standard deviation of speed
%
%    varargin - (Optional) name value pair to be set for eye movement. See
%               emSet for supported parameter list.
%
% Outputs:
%   em        - eye movement structure
%
% Optional key/value pairs:
%    **Needs to be filled out*
%
% Notes:
%    1) For all eye-movements, we assume that the acceleration time is
%       negligible
%    2) Drift and tremor can be superimposed. However, during the
%       microsaccade, both drift and tremor are suppressed
%    3) We assume that drift is a 2D brownian motion. This is
%       reasonable when micro-saccade is present. However, when
%       micro-saccade is suppressed, drift will somehow towards the
%       fixation point. The actual pattern is unknown and this
%       return-to-origin feature is not implemented here.
%    4) For micro-saccade, speed and amplitude follows main-squence as
%       saccade and here, we just use speed
%    5) We don't add a field 'duration' to microsaccade and in computation, 
%       the duration of microsaccade will be estimated by the speed and
%       distance between current position to the fixation point
%
% References:
%    1) Susana Martinez-Conde et. al, The role of fixational eye movements
%       in visual perception, Nature reviews | neuroscience, Vol. 5, 2004, 
%       page 229~240
%    2) Susana Martinez-Conde et. al, Microsaccades: a neurophysiological
%       analysis, Trends in Neurosciences, Volume 32, Issue 9, September
%       2009, Pages 463~475
%
% See Also:
%    emGenSequence, emSet, emGet
%

% History:
%    xx/xx/14  HJ   Copyright ISETBIO TEAM 2014
%    05/02/18  jnm  Formatting

% Examples:
%{
	em = emCreate;
%}

% Fill in default values for missing fields in params
if notDefined('params'), params = []; end

% set params to default values
% set general fields
p.name = 'em structure';
p.type = 'eye movement';
p.emFlag = ones(3, 1); % emType - no eye movement
p.sampTime = 0.001; % sample time interval - 1 ms

% set fields for tremor
p.tremor.interval = 0.012;    % Tremor mean frequency - 83 Hz
p.tremor.intervalSD = 0.001;  % std of tremor frequency - 60~100 Hz
p.tremor.amplitude = 0.0073;  % Tremor amplitude -  18 arcsec

% set fields for drift
%   There's a big difference for drift speed between literatures, we just
%   pick a reasonable value among them
p.drift.speed = 3 / 60 * pi / 180;    % drift speed - drift speed, rad/sec
p.drift.speedSD = 1 / 60 * pi / 180;  % std of drift speed

% set fields for micro-saccades
p.msaccade.interval = 0.6;          % micro-saccade interval - 0.6 secs
p.msaccade.intervalSD = 0.3;        % std for micro-saccade interval
p.msaccade.dirSD = 5 * pi / 180;    % std for direction
p.msaccade.speed = 15 * pi / 180;   % micro saccade speed - 15 deg/s
p.msaccade.speedSD = 5 * pi / 180;  % std for micro saccade speed

% merge params with default values
if ~isfield(params, 'msaccade') && isfield(params, 'microsaccade')
    params.msaccade = params.microsaccade;
    params = rmfield(params, 'microsaccade');
end

% Cool, but I am worried about field name incompatibilities
% Maybe we should do this as a for loop with emSet(...)
% And also, I couldn't find the function so I asked HJ where he got it.
% em = setstructfields(p, params);
% process values in params
em = overwriteStruct(p, params);

% process values in varargin
for ii = 1:2:length(varargin)
    em = emSet(em, varargin{ii}, varargin{ii + 1});
end

% some checks for params
assert(numel(em.emFlag) == 3, 'emType should be 3x1 logical vector');

end

%% Set structfields in our own form
function params = overwriteStruct(params, newP)
% Replace field values in params with the values in newP
%
% Syntax:
%   params = overwriteStruct(params, newP)
%
% Description:
%    There is another implementation for this function in recent Matlab
%    release in signal processing toolbox
%
% Inputs:
%    params - Struct. Starting structure
%    newP   - Struct. Structure containing values to replace param's vals.
%
% Outputs:
%    params - Struct. The modified params structure.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/15  HJ   ISETBIO TEAM, 2015
%    05/02/18  jnm  Formatting

% check inputs
if ~isstruct(params) || ~isstruct(newP), return; end

% Loop over all the fields of params
fields = fieldnames(params);
for ii = 1:length(fields)
    if isfield(newP, fields{ii})
        v = newP.(fields{ii});
        if isstruct(params.(fields{ii})) && isstruct(v)
            params.(fields{ii}) = overwriteStruct(params.(fields{ii}), v);
        else
            params.(fields{ii}) = v;
        end
    end
end

end