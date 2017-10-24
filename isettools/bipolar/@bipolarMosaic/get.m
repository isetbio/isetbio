function val = get(obj, varargin)
%  Gets isetbio bipolar object parameters
% 
%    val = bipolar.get(parameter)
% 
% Parameters
%     cellLocation;               % location of bipolar RF center
%     patchSize;                  % size of retinal patch from sensor
%     timeStep;                   % time step of simulation from sensor
%     sRFcenter;                  % spatial RF of the center on the receptor grid
%     sRFsurround;                % spatial RF of the surround on the receptor grid
%     temporalDifferentiator;     % differentiator function
%     responseCenter;             % Store the linear response of the center after convolution
%     responseSurround;           % Store the linear response of the surround after convolution
%     response;
%     threshold;
%     rgbdata;
%
% 5/2016 JRG (c) isetbio team
%%
p = inputParser;
p.CaseSensitive = false; 
p.FunctionName = mfilename;
p.KeepUnmatched = true;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFields = {...
    'cellLocation',...
    'patchSize',...
    'timeStep',...
    'sRFcenter',...
    'sRFsurround',...
    'spatialrf',...
    'temporalDifferentiator',...
    'threshold',...
    'rgbdata',...
    'responseCenter','bipolarresponsecenter',...
    'responseSurround','bipolarresponsesurround',...
    'response','bipolarresponse'...
    'duration' ...
    };

p.addRequired('what',@(x) any(validatestring(ieParamFormat(x),allowableFields)));

% Parse and put results into structure p.
p.parse(varargin{:}); params = p.Results;

switch ieParamFormat(params.what)  % Lower case and remove spaces

    case {'celllocation'}        
        val = obj.cellLocation;
                
    case{'patchsize'}
        val = obj.patchSize;
        
    case{'timestep'}
        val = obj.timeStep;
        
    case {'duration'}
        % In seconds
        val = obj.timeStep*size(obj.responseCenter,3);
        
    case{'rgbdata'}
        val = obj.rgbData;
        
    case{'srfcenter'}
        val = obj.sRFcenter;
        
    case{'srfsurround'}
        val = obj.sRFsurround;
        
    case {'spatialrf'}
        val = obj.sRFcenter - obj.sRFsurround;
        
    case{'temporaldifferentiator'}
        val = obj.temporalDifferentiator;
        
    case{'threshold'}
        val = obj.threshold;
        
    case{'responsecenter','bipolarresponsecenter'}
        val = obj.responseCenter;
        
    case{'responsesurround','bipolarresponsesurround'}
        val = obj.responseSurround;
        
    case{'response','bipolarresponse'}
        val = obj.responseCenter - obj.responseSurround;
end

