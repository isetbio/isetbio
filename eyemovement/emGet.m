function val = emGet(em, params, varargin)
% Get eye movement properties
%
% Syntax:
%   val = emGet(em, params, [varargin])
%
% Description:
%    The model for generating an eye movement sequence is defined in
%    @coneMosaic.emGenSequence.
%
%    Retrieve the parameter specified in params.
%
% Inputs:
%    em       - Struct. Eye movement structure, created by emCreate
%    params   - String. Parameter name to get, spaces and case insensitive.
%               Options include the following:
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
%       {'drift speed'}          - Numeric. Drift speed. Default is
%                                  rad/sec, but cones/sample time is also
%                                  an acceptable unit.
%       {'drift speed SD'}       - Numeric. Std Dev. of drift speed
%
%       {'msaccade interval'}    - Numeric. Time between micro-saccade
%       {'msaccade interval SD'} - Numeric. Std Dev. of m-saccade interval
%       {'msaccade dir SD'}      - Numeric. Direction variation in radians.
%       {'msaccade speed'}       - Numeric. Micro-saccade speed
%       {'msaccade speed SD'}    - Numeric. Micro-saccade speed std dev.
%    varargin - (Optional) Additional information that may be required to
%               retrieve the desired parameter. One possible addition would
%               be units for some parameters.
%
% Outputs:
%    val      - value for parameter, if not found, return empty
%
% Notes:
%    * The default drift speed in terms of cones per second is
%        emGet(emCreate, 'drift speed', 'cones/sample') * 1e+3  ~10
%      This is slower than Rucci/Victor argue for in their 2015 TICS paper,
%      where they say for natural vision the speed is 50'/sec which is
%      about 100 cones/sec. The value here is consistent with the value
%      reported for experienced psychophysical observers who fixate well
%      (also according to Rucci and Victor, citing Ditchburn but also
%      Susana Martinez Condes paper.
%
% See Also:
%    emSet, emCreate, @coneMosaic.emGenSequence
%

% History:
%    xx/xx/14  HJ/BW  (c) ISETBIO Team, 2014
%    04/19/18  jnm    Formatting & minor fixes.

%  Examples:
%{
    em = emCreate;
    fs = emGet(em, 'frequency');
    params.f = 0.017;
    params.w = 1.5e-6;
    amp = emGet(em, 'tremor amplitude', 'cones', params);
%}

%% Check inputs
if notDefined('em'), error('eye movement structure required'); end
if notDefined('params'), error('Parameter name required'); end

%% Get property value
val = [];
params = ieParamFormat(params);  % Lower case and remove spaces
switch params
    % Basic information about eye movement
    case {'name'}
        % name = emGet(em, 'name');
        if isfield(em, 'name'), val = em.name; end
    case {'type'}
        % type = emGet(em, 'type');
        if isfield(em, 'type'), val = em.type; end
    case {'emflag', 'flag'}
        % emFlag = emGet(em, 'em flag');
        if isfield(em, 'emFlag'), val = em.emFlag; end
    case {'sampletime', 'samptime'}
        % sampTime = emGet(em, 'sample time', unit);
        if isfield(em, 'sampTime'), val = em.sampTime; end
        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1});
        end
    case {'frequency', 'samplefrequency'}
        % fs = emGet(em, 'frequency');
        % Frequency in HZ
        val = 1 / emGet(em, 'sample time');

    % Eye movement substructures - tremor, drift, micro-saccade
    case {'tremor'}
        % tremor = emGet(em, 'tremor');
        if isfield(em, 'tremor'), val = em.tremor; end
    case {'drift'}
        % drift = emGet(em, 'drift');
        if isfield(em, 'drift'), val = em.drift; end
    case {'microsaccade', 'msaccade'}
        % msaccade = emGet(em, 'msaccade');
        if isfield(em, 'msaccade'), val = em.msaccade; end

    % Tremor parameters
    %    The tremor sample times are spaced on average by tremor interval
    %    with a standard devation of tremoral interval sd. The standard
    %    deviation of the tremor size is stored in tremoral interval std
    %    dev. The value of this standard deviation is stored per second, so
    %    the standard deviaton for a particular tremor interval is
    %    (interval * std dev).
    case {'tremorinterval'}
        % interval = emGet(em, 'tremor interval');
        if checkfields(em, 'tremor', 'interval')
            val = em.tremor.interval;
        end
        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1});
        end
    case {'tremorintervalsd'}
        % intervalSD = emGet(em, 'tremor interval SD');
        if checkfields(em, 'tremor', 'intervalSD')
            val = em.tremor.intervalSD;
        end
        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1});
        end
    case {'tremoramplitude'}
        % amp = emGet(em, 'tremor amplitude');
        if checkfields(em, 'tremor', 'amplitude')
            % Default unit is radians
            val = em.tremor.amplitude;
        end
        if length(varargin) > 1
            val = val * emUnitScaleFactor(em, varargin{1}, varargin{2});
        elseif ~isempty(varargin)
            val = val * emUnitScaleFactor(em, varargin{1});
        end

    % Drift parameters
    case {'driftspeed'}
        % speed = emGet(em, 'drift speed');
        if checkfields(em, 'drift', 'speed'), val = em.drift.speed; end
        if length(varargin) > 1
            val = val * emUnitScaleFactor(em, varargin{1}, varargin{2});
        elseif ~isempty(varargin)
            val = val * emUnitScaleFactor(em, varargin{1});
        end
    case {'driftspeedsd'}
        % speedSD = emGet(em, 'drift speed SD');
        if checkfields(em, 'drift', 'speedSD')
            val = em.drift.speedSD;
        end
        if length(varargin) > 1
            val = val * emUnitScaleFactor(em, varargin{1}, varargin{2});
        elseif ~isempty(varargin)
            val = val * emUnitScaleFactor(em, varargin{1});
        end

    % Micro-saccade parameters
    case {'msaccadeinterval', 'microsaccadeinterval'}
        % interval = emGet(em, 'msaccade interval');
        if checkfields(em, 'msaccade', 'interval')
            val = em.msaccade.interval;
        end
        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1});
        end
    case {'msaccadeintervalsd', 'microsaccadeintervalsd'}
        % intervalSD = emGet(em, 'msaccade interval SD');
        if checkfields(em, 'msaccade', 'intervalSD')
            val = em.msaccade.intervalSD;
        end
        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1});
        end
    case {'msaccadedirsd', 'microsaccadedirsd'}
        % dirSD = emGet(em, 'msaccade dir SD');
        if checkfields(em, 'msaccade', 'dirSD')
            val = em.msaccade.dirSD;
        end
        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1});
        end
    case {'msaccadespeed', 'microsaccadespeed'}
        % speed = emGet(em, 'msaccade speed'); 
        if checkfields(em, 'msaccade', 'speed')
            val = em.msaccade.speed;
        end
        if length(varargin) > 1
            val = val * emUnitScaleFactor(em, varargin{1}, varargin{2});
        elseif ~isempty(varargin)
            val = val * emUnitScaleFactor(em, varargin{1});
        end
    case {'msaccadespeedsd', 'microsaccadespeedsd'}
        % speedSD = emGet(em, 'msaccade speed SD'); 
        if checkfields(em, 'msaccade', 'speedSD')
            val = em.msaccade.speedSD;
        end
        if length(varargin) > 1
            val = val * emUnitScaleFactor(em, varargin{1}, varargin{2});
        elseif ~isempty(varargin)
            val = val * emUnitScaleFactor(em, varargin{1});
        end

    % Catch-all for unknown/invalid parameters.
    otherwise
        error('Unknown parameter encountered');
end

end

%% Aux-function: Get unit scale factor for eye movement
function val = emUnitScaleFactor(em, unitName, params)
% Converts from standard units to other scales
%
% Syntax:
%   val = emUnitScaleFactor(em, unitName, params)
%
% Description:
%    Converts from standard units to other scales.
%
% Inputs:
%    em       - Struct. An eye movement structure.
%    unitName - String. The parameter's unit name.
%    params   - Struct. A structure containing frequency and wavelength.
%
% Outputs:
%    val      - Numeric. The units in the desired scale.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    ieUnitScaleFactor
%

% History:
%    xx/xx/14  HJ   ISETBIO Team, 2014
%    04/19/18  jnm  Formatting

% Init
if notDefined('em'), em = emCreate; end
if notDefined('unitName'), error('Unit name required'); end
if notDefined('params'), params = []; end
unitName = ieParamFormat(unitName);

% Get unit scale factor
switch unitName
    case {'sampletimes', 'samplingtime'} % for time
        % From secs to number of sampling times
        val = emGet(em, 'frequency');
    case {'cones/sample', 'conespersample'} % for speed
        % From rad/sec to number of cones per sample time
        % params.f contains focal length in meters (default 0.017)
        % params.w contains cone width (default 1.5e-6)
        if notDefined('params.f'), f = 0.017; else, f = params.f; end
        if notDefined('params.w'), w = 1.5e-6; else, w = params.w; end
        mperdeg  = tand(1) * f;
        sampTime = emGet(em, 'sample time');
        val = 180/pi * sampTime * mperdeg / w;
    case {'cones', 'samples'} % distance, amplitude
        % From rad to number of cones
        % params.f contains focal length in meters (default 0.017)
        % params.w contains cone width (default 1.5e-6)
        if notDefined('params.f'), f = 0.017; else, f = params.f; end
        if notDefined('params.w'), w = 1.5e-6; else, w = params.w; end
        mperdeg = tand(1) * f;
        val = 180/pi * mperdeg / w;
    otherwise
        val = ieUnitScaleFactor(unitName);
end

end