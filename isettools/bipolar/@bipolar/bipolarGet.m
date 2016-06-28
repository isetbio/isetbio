function val = bipolarGet(obj, varargin)
%  Gets isetbio bipolar object parameters.
% 
% The bipolar object allows the simulated cone responses to be passed on to
% the inner retina object and to approxiately maintain its impulse
% response. This will allow us to run the nonlinear biophysical cone outer
% segment model and pass its results on to the bipolar stage and then RGCs.
% 
%     cellLocation;                    % location of bipolar RF center
%     patchSize;                       % size of retinal patch from sensor
%     timeStep;                        % time step of simulation from sensor
%     sRFcenter;                       % spatial RF of the center on the receptor grid
%     sRFsurround;                     % spatial RF of the surround on the receptor grid
%     temporalDifferentiator;          % differentiator function
%     responseCenter;                  % Store the linear response of the center after convolution
%     responseSurround;                % Store the linear response of the surround after convolution
% 
% 
% 5/2016 JRG (c) isetbio team
%%
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
    'temporalDifferentiator',...
    'threshold',...
    'responseCenter','bipolarresponsecenter',...
    'responseSurround','bipolarresponsesurround',...
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
        
    case{'srfcenter'}
        val = obj.sRFcenter;
        
    case{'srfsurround'}
        val = obj.sRFsurround;
        
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

