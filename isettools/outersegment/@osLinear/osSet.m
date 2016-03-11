function obj = osSet(obj, varargin)
% Sets isetbio outersegment object parameters
% 
%  One parameter can be set by each call
%
% Parameters:
%   noiseFlag   -  sets current as noise-free ('0') or noisy ('1')
%   sConeFilter - the linear filter for S-cone temporal response
%   mConeFilter - the linear filter for M-cone temporal response
%   lConeFilter - the linear filter for L-cone temporal response
%   conespacing  -
%   conesampling -
%   conecurrentsignal - 
%
% Example:
%   
%   adaptedOS = osSet(adaptedOS, 'noiseFlag', 0);
% 
% 8/2015 JRG NC DHB

%% Check for the number of arguments and create parser object.
% 
% Check key names with a case-insensitive string, errors in this code are
% attributed to this function and not the parser object.
narginchk(0, Inf);
p = inputParser; 
p.CaseSensitive = false; 
p.FunctionName = mfilename;

% Force to lower case and eliminate spaces
for ii=1:2:length(varargin), varargin{ii} = ieParamFormat(varargin{ii}); end

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFieldsToSet = {...
    'noiseflag',...
    'sconefilter',...
    'mconefilter',...
    'lconefilter',...    
    'conespacing',...
    'conesampling',...
    'conecurrentsignal'};
p.addRequired('what',@(x) any(validatestring(x,allowableFieldsToSet)));
p.addRequired('value');

% Parse and put results into structure p.
p.parse(varargin{:}); 
params = p.Results;

%% Set key-value pairs.
switch lower(params.what)

    case{'noiseflag'}
        % noise flag = 0
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
        
    case{'conespacing'}
        % Spatial sample spacing
        obj.coneSpacing = params.value;
        
    case{'conesampling'}
        % Temporal sample spacing
        obj.coneSampling = params.value;
        
    case{'conecurrentsignal'}
        obj.coneCurrentSignal = params.value; 
        
end

end