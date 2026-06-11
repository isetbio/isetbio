function obj = osSet(obj, varargin)
% Sets isetbio outersegment object parameters.
%
% Syntax:
%   obj = osSet(obj, varargin)
%
% Description:
%    Set the isetbio outersegment object parameter to the provided value.
%
%    Examples are contained in the code. To access, type 'edit osSet.m'
%    into the Command Window.
%
% Inputs:
%    obj      - The outersegment object
%    varargin - Additional information required to perform the specified
%               operations. Parameters to set include the following:
%        {'sConeFilter'} - the linear filter for S-cone temporal response
%        {'mConeFilter'} - the linear filter for M-cone temporal response
%        {'lConeFilter'} - the linear filter for L-cone temporal response
%        {'noiseFlag'}   -  sets current as noise-free ('0') or noisy ('1')
%        {'patchSize'}   - cone current as a function of time
%        {'timeStep'}    - noisy cone current signal
%        {'size'}        - array size of photon rate
%        {'coneCurrentSignal'}
%                        - cone current as a function of time
%
% Outputs:
%    obj      - The modified outersegment object
%
% Optional key/value pairs:
%    None.
%
% Notes:
%

% History:
%    08/xx/15  JRG NC DHB  Created
%    02/13/18  jnm         Formatting

% Examples:
%{
    adaptedOS = osCreate;
	adaptedOS = osSet(adaptedOS, 'noiseFlag', 'none');
    adaptedOS = osSet(adaptedOS,  'time step', 0.001);
%}

narginchk(0, Inf);
p = inputParser;
p.CaseSensitive = false;
p.FunctionName = mfilename;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFieldsToSet = {...
    'noiseflag', ...
    'sconefilter', ...
    'mconefilter', ...
    'lconefilter', ...
    'patchsize', ...
    'timestep', ...
    'size', ...
    'conecurrentsignal'};
p.addRequired('what', @(x) any(validatestring(ieParamFormat(x), ...
    allowableFieldsToSet)));
p.addRequired('value');

% Parse and put results into structure p.
p.parse(varargin{:});
params = p.Results;

%% Set all the key-value pairs.

    switch ieParamFormat(params.what)  % Lower case and remove spaces

        case {'noiseflag'}
            if ischar(params.value) && (ismember(lower(params.value), ...
                    outerSegment.validNoiseFlags))
                obj.noiseFlag = params.value;
            else
                error('%s is an illegal value for os.noiseFlag', ...
                    params.value);
            end

        case {'sconefilter'}
            % Temporal impulse response
            obj.lmsConeFilter(:, 3) = value;

        case {'mconefilter'}
            % Temporal impulse response
            obj.lmsConeFilter(:, 2) = value;

        case {'lconefilter'}
            % Temporal impulse response
            obj.lmsConeFilter(:, 1) = value;

        case{'patchsize'}
            obj.patchSize = params.value;

        case{'timestep'}
            obj.timeStep = params.value;

        case{'conecurrentsignal'}
            obj.coneCurrentSignal = params.value;
    end
    