function em = emSet(em, params, val)
% Assign values to eye movement properties
%
% Syntax:
%   em = emSet(em, params, val)
%
% Description:
%    Assign values to eye movement properties.
%
%  Inputs:
%    em       - Struct. eye movement structure, created by emCreate
%    params   - String. parameter name to get, space and case insensitive.
%               Value type is indicated in the description. Options include
%               the following:
%       {'name'}                 - String. User defined eye movement name.
%       {'em flag'}              - Vector. 3x1 vector, indicating whether
%                                  or not to include tremor, drift and
%                                  micro-saccade.
%       {'sample time'}          - Numeric. The number of time units per
%                                  sample. 'sec', or 'ms'. Default 'sec'.
%       {'frequency'}            - Numeric. Sampling frequency (Hz)
%       {'tremor'}               - Struct. Tremor structure.
%       {'drift'}                - Struct. Drift structure.
%       {'micro saccade'}        - Struct. Micro-saccade structure.
%
%       {'tremor interval'}      - Numeric. Time between tremors in secs.
%       {'tremor interval SD'}   - Numeric. Std Dev. of tremor interval.
%       {'tremor amplitude'}     - Numeric. Amplitude of tremor in radians.
%
%       {'drift speed'}          - Numeric. Drift speed (rad/sec)
%       {'drift speed SD'}       - Numeric. Std Dev. of drift speed
%
%       {'msaccade interval'}    - Numeric. Time between micro-saccade
%       {'msaccade interval SD'} - Numeric. Std Dev. of m-saccade interval
%       {'msaccade dir SD'}      - Numeric. Direction variation in radians.
%       {'msaccade speed'}       - Numeric. Micro-saccade speed
%       {'msaccade speed SD'}    - Numeric. Micro-saccade speed std dev.
%    val      - value to be set to em.(params). Type will depend on the
%               parameter, and can be found above with the parameter.
%
%  Outputs:
%    em       - Struct. The modified eye movement structure.
%
%  See also:
%    emSet, emCreate
%

% History:
%    xx/xx/14  HJ/BW  (c) ISETBIO Team, 2014
%    04/19/18  jnm    Formatting, fixed example. Added listed examples that
%                     were not actually present in switch case (frequency)

%  Examples:
%{
    em = emCreate;
    fs = emGet(em, 'frequency')

    em = emSet(em, 'frequency', fs*2)
    fs = emGet(em, 'frequency')
%}

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
    case {'emflag', 'flag'}
        assert(numel(val) == 3, 'emFlag should be 3x1 vector');
        em.emFlag = val(:);
    case {'sampletime'}
        em.sampTime = val;
    case {'frequency', 'samplefrequency'}
        em.sampTime = 1/val;

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