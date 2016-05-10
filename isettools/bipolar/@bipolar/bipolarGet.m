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
    'sRFcenter',...
    'sRFsurround',...
    'tIR',...
    'temporalDifferentiator',...
    'threshold',...
    'responseCenter','bipolarresponsecenter',...
    'responseSurround','bipolarresponsesurround'...
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
        
    case{'srfcenter'}
        val = obj.sRFcenter;
        
    case{'srfsurround'}
        val = obj.sRFsurround;
        
    case{'tIR'}
        val = obj.tIR;
        
    case{'temporaldifferentiator'}
        val = obj.temporalDifferentiator;
        
    case{'threshold'}
        val = obj.threshold;
        
    case{'responsecenter','bipolarresponsecenter'}
        val = obj.responseCenter;
        
    case{'responsesurround','bipolarresponsesurround'}
        val = obj.responseCenter;
end

