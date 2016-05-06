function val = bipolarGet(obj, varargin)
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
    'sRF',...
    'tIR',...
    'threshold',...
    'response','bipolarresponse'...
    };
p.addRequired('what',@(x) any(validatestring(ieParamFormat(x),allowableFieldsToSet)));

% Parse and put results into structure p.
p.parse(varargin{:}); params = p.Results;

switch ieParamFormat(params.what);  % Lower case and remove spaces

    case {'celllocation'}        
        val = obj.cellLocation;
                
    case{'patchsize'}
        val = obj.patchSize;
        
    case{'timestep'}
        val = obj.timeStep;
        
    case{'rgbdata'}
        val = obj.rgbData;
        
    case{'srf'}
        val = obj.sRF;
        
    case{'tIR'}
        val = obj.tIR;
        
    case{'threshold'}
        val = obj.threshold;
        
    case{'response','bipolarresponse'}
        val = obj.response;
        
end

