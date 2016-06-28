function obj = osSet(obj, varargin)
% Sets isetbio outersegment object parameters.
%
% Parameters:
% 
%       {'sConeFilter'} - the linear filter for S-cone temporal response
%       {'mConeFilter'} - the linear filter for M-cone temporal response
%       {'lConeFilter'} - the linear filter for L-cone temporal response
%       {'noiseFlag'} -  sets current as noise-free ('0') or noisy ('1')
%       {'patchSize'} - cone current as a function of time
%       {'timeStep'} - noisy cone current signal
%       {'size'} - array size of photon rate
%       {'coneCurrentSignal'} - cone current as a function of time
% 
% Example:
%
%   adaptedOS = osSet(adaptedOS, 'noiseFlag', 0, 'time step',0.001);
%
% 8/2015 JRG NC DHB

narginchk(0, Inf);
p = inputParser; p.CaseSensitive = false; p.FunctionName = mfilename;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFieldsToSet = {...
    'noiseflag',...
    'sconefilter',...
    'mconefilter',...
    'lconefilter',...
    'patchsize',...
    'timestep',...
    'size',...
    'conecurrentsignal'};
p.addRequired('what',@(x) any(validatestring(ieParamFormat(x),allowableFieldsToSet)));
p.addRequired('value');

% Parse and put results into structure p.
p.parse(varargin{:}); params = p.Results;


%% Set all the key-value pairs.

    switch ieParamFormat(params.what);  % Lower case and remove spaces
        
        case {'noiseflag'}
            obj.noiseFlag = params.value;
            
        case {'sconefilter'}
            % Temporal impulse response
            obj.sConeFilter = params.value;
            
        case {'mconefilter'}
            % Temporal impulse response
            obj.mConeFilter = params.value;
            
        case {'lconefilter'}
            % Temporal impulse response
            obj.lConeFilter = params.value;
            
        case{'patchsize'}
            obj.patchSize = params.value;
            
        case{'timestep'}
            obj.timeStep = params.value;
            
        case{'conecurrentsignal'}
            obj.coneCurrentSignal = params.value;
    end
    