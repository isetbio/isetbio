function em = emSet(em, params, val)
%% function emSet(em, params, val)
%    Get properties from cones
%
%  Inputs:
%    em       - eye movement structure, created by emCreate
%    params   - parameter name to get, spaces and case insensitive
%    val      - value to be set to em.(params)
%
%  Outputs:
%    em       - eye movement structure with params set
%
%  Supported params:
%    {'name'}                 - user defined name of eye movement
%    {'em flag'}              - 3x1 vector, indicating whether to include
%                               tremor, drift and micro-saccade
%    {'sample time'}          - secs per sample ([sec], ms)
%    {'frequency'}            - sampling frequency (Hz)
%
%    {'tremor'}               - structure for tremor
%    {'drift'}                - structure for drift
%    {'micro saccade'}        - structure for micro-saccade
%
%    {'tremor interval'}      - time between tremors in secs
%    {'tremor interval SD'}   - standard deviation of tremor interval
%    {'tremor amplitude'}     - amplitude of tremor in rads
%
%    {'drift speed'}          - drift speed (rad/sec)
%    {'drift speed SD'}       - standard deviation of drift speed
%
%    {'msaccade interval'}    - time between micro-saccade
%    {'msaccade interval SD'} - standard deviation of m-saccade interval
%    {'msaccade dir SD'}      - direction variation
%    {'msaccade speed'}       - micro-saccade speed
%    {'msaccade speed SD'}    - standard deviation for micro-saccade speed
%
%  Example:
%    em = emCreate;
%    fs = emGet(em, 'frequency');
%    
%    amp = emGet(em, 'tremor amplitude');
%
%  See also:
%    emSet, emCreate
%
%
%  HJ/BW (c) ISETBIO Team, 2014

%% Check inputs
if notDefined('em'), error('eye movement structure required'); end
if notDefined('params'), error('Parameter name required'); end
if notDefined('val'), error('Value required'); end

%% Get property value
params = ieParamFormat(params);  % Lower case and remove spaces
switch params
    % Basic information about eye movement
    case {'name'}
        if ischar(val), em.name = val; end
    case {'emflag'}
        assert(numel(val) == 3, 'emFlag should be 3x1 vector');
        em.emFlag = val(:);
    case {'sampletime'}
        em.sampTime = val;
        
    % Eye movement substructures - tremor, drift, micro-saccade
    case {'tremor'}
        assert(isstruct(val), 'val should be a structure');
        em.tremor = val;
    case {'drift'}
        assert(isstruct(val), 'val should be a structure');
        em.drift = val;
    case {'microsaccade', 'msaccade'}
        assert(isstruct(val), 'val should be a structure');
        em.msaccade = val;
        
    % tremor parameters
    case {'tremorinterval'}
        % em = emSet(em, 'tremor interval', val);
        em.tremor.interval = val;
    case {'tremorintervalsd'}
        em.tremor.intervalSD = val;
    case {'tremoramplitude'}
        em.tremor.amplitude = val;
    
    % drift parameters
    case {'driftspeed'}
        em.drift.speed = val;
    case {'driftspeedsd'}
        em.drift.speedSD = val;
      
    % micro-saccade parameters
    case {'msaccadeinterval', 'microsaccadeinterval'}
        em.msaccade.interval = val;
    case {'msaccadeintervalsd', 'microsaccadeintervalsd'}
        em.msaccade.intervalSD = val;
    case {'msaccadedirsd', 'microsaccadedirsd'}
        em.msaccade.dirSD = val;
    case {'msaccadespeed', 'microsaccadespeed'}
        em.msaccade.speed = val;
    case {'msaccadespeedsd', 'microsaccadespeedsd'}
        em.msaccade.speedSD = val;
    otherwise
        error('Unknown parameter encountered');
end

end