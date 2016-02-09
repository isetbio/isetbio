function obj = osSet(obj, varargin)
% Sets isetbio outersegment object parameters using the input parser structure.
% 
% Parameters:
%       {'noiseFlag'} -  sets current as noise-free ('0') or noisy ('1')
%       {'sConeFilter'} - the linear filter for S-cone temporal response
%       {'mConeFilter'} - the linear filter for M-cone temporal response
%       {'lConeFilter'} - the linear filter for L-cone temporal response
% 
% adaptedOS = osSet(adaptedOS, 'noiseFlag', 0);
% 
% 8/2015 JRG NC DHB

% Check for the number of arguments and create parser object.
% Parse key-value pairs.
% 
% Check key names with a case-insensitive string, errors in this code are
% attributed to this function and not the parser object.
narginchk(0, Inf);
p = inputParser; p.CaseSensitive = false; p.FunctionName = mfilename;

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

% Define what units are allowable.
% allowableUnitStrings = {'a', 'ma', 'ua', 'na', 'pa'}; % amps to picoamps

% Set up key value pairs.
% Defaults units:
% p.addParameter('units','pa',@(x) any(validatestring(x,allowableUnitStrings)));
% p.addParameter('sconefilter',0,@isnumeric);
% p.addParameter('mconefilter',0,@isnumeric);
% p.addParameter('lconefilter',0,@isnumeric);

% Parse and put results into structure p.
p.parse(varargin{:}); params = p.Results;

% % Old error check on input.
% if ~exist('params','var') || isempty(params)
%     error('Parameter field required.');
% end
% if ~exist('val','var'),   error('Value field required.'); end;

% Set key-value pairs.
switch lower(params.what)

    case{'noiseflag'}
        obj.noiseFlag = params.value;
        
    case {'sconefilter'}
        obj.sConeFilter = params.value;
        
    case {'mconefilter'}
        obj.mConeFilter = params.value;
        
    case {'lconefilter'}
        obj.lConeFilter = params.value;
        
    case{'conespacing'}
        obj.coneSpacing = params.value;
        
    case{'conesampling'}
        obj.coneSampling = params.value;
        
    case{'conecurrentsignal'}
        obj.coneCurrentSignal = params.value; 
        
end

