function val = emGet(em, params, varargin)
%% function emGet(em, params, [varargin])
%    Get properties from cones
%
%  Inputs:
%    em       - eye movement structure, created by emCreate
%    params   - parameter name to get, spaces and case insensitive
%    varargin - possible units for some parameters
%
%  Outputs:
%    val      - value for parameter, if not found, return empty
%
%  Supported params:
%    {'name'}                 - user defined name of eye movement
%    {'type'}                 - should be 'eye movement'
%    {'em flag'}              - 3x1 vector, indicating whether to include
%                               tremor, drift and micro-saccade
%    {'sample time'}          - secs per sample ([sec], ms)
%    {'frequency'}            - sampling frequency (Hz)
%
%    {'tremor'}               - structure for tremor
%    {'drift'}                - structure for drift
%    {'micro saccade'}        - structure for micro-saccade
%
%    {'tremor interval'}      - time between tremors ([sec], ms, samples)
%    {'tremor interval SD'}   - standard deviation of tremor interval
%    {'tremor amplitude'}     - amplitude of tremor ([rad], deg, cones)
%
%    {'drift speed'}          - drift speed ([rad/s], cones/sample time)
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
%    params.f = 0.017; params.w = 1.5e-6;
%    amp = emGet(em, 'tremor amplitude', 'cones', params);
%
%  See also:
%    emSet, emCreate
%
%
%  HJ/BW (c) ISETBIO Team, 2014

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
    case {'emflag'}
        % emFlag = emGet(em, 'em flag');
        if isfield(em, 'emFlag'), val = em.emFlag; end
    case {'sampletime'}
        % sampTime = emGet(em, 'sample time',unit);
        % 
        if isfield(em, 'sampTime'), val = em.sampTime; end
        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1});
        end
    case {'frequency', 'samplefrequency'}
        % fs = emGet(em, 'frequency');
        % Frequency in HZ
        val = 1/emGet(em,'sample time');

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
        
    % tremor parameters
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
            val = em.tremor.amplitude;
        end
        if length(varargin) > 1
            val = val * emUnitScaleFactor(em, varargin{1}, varargin{2});
        elseif ~isempty(varargin)
            val = val * emUnitScaleFactor(em, varargin{1});
        end
    
    % drift parameters
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
      
    % micro-saccade parameters
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
    
    otherwise
        error('Unknown parameter encountered');
end

end

%% Aux-function: Get unit scale factor for eye movement
function val = emUnitScaleFactor(em, unitName, params)
% Converts from standard units to other scales
%
% See also:
%   ieUnitScaleFactor
%
% (HJ) ISETBIO Team, 2014

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
        if notDefined('params.f'), f = 0.017;  else f = params.f; end
        if notDefined('params.w'), w = 1.5e-6; else w = params.w; end
        mperdeg  = tand(1) * f;
        sampTime = emGet(em, 'sample time');
        val = 180/pi * sampTime * mperdeg / w;
    case {'cones', 'samples'} % distance, amplitude
        % From rad to number of cones
        % params.f contains focal length in meters (default 0.017)
        % params.w contains cone width (default 1.5e-6)
        if notDefined('params.f'), f = 0.017;  else f = params.f; end
        if notDefined('params.w'), w = 1.5e-6; else w = params.w; end
        mperdeg = tand(1) * f;
        val = 180/pi * mperdeg / w;
    otherwise
        val = ieUnitScaleFactor(unitName);
end

end