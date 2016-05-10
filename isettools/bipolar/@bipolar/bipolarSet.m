function obj = bipolarSet(obj, varargin)
%  Gets isetbio bipolar object parameters.
% 
% Parameters:
%       {''} -  
% 
% 5/2016 JRG 

narginchk(0, Inf);
p = inputParser; p.CaseSensitive = false; p.FunctionName = mfilename;
p.KeepUnmatched = true;
% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFieldsToSet = {...
    'cellLocation',...
    'patchSize',...
    'timeStep',...
    'sRFcenter',...
    'sRFsurround',...
    'tIR',...
    'temporalDifferentiator',...
    'threshold',...
    'response'...
    };
p.addRequired('what',@(x) any(validatestring(ieParamFormat(x),allowableFieldsToSet)));
p.addRequired('value');

% Parse and put results into structure p.
p.parse(varargin{:}); params = p.Results;

switch ieParamFormat(params.what);  % Lower case and remove spaces

    case {'celllocation'}        
        obj.cellLocation = params.value;
                
    case{'patchsize'}
        obj.patchSize = params.value;
        
    case{'timestep'}
        obj.timeStep = params.value;
        
    case{'rgbdata'}
        obj.rgbData = params.value;
        
    case{'srfcenter'}
        obj.sRFcenter = params.value;
        
    case{'srfsurround'}
        obj.sRFsurround = params.value;
        
    case{'tIR'}
        obj.tIR = params.value;
        
    case{'threshold'}
        obj.threshold = params.value;
        
    case{'response'}
        obj.response = params.value;
        
end

