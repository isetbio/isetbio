function obj = bipolarSet(obj, varargin)
% Assigns value to isetbio bipolar object parameters.
%
% Syntax:
%
%   bp.bipolarSet('PARAM1', val1, 'PARAM2', val2);
%
% Description:
%    Sets isetbio bipolar object parameters.
% 
%    The bipolar object allows the simulated cone responses to be passed on
%    to the inner retina object and to approxiately maintain its impulse
%    response. This will allow us to run the nonlinear biophysical cone
%    outer segment model and pass its results on to the bipolar stage and
%    then RGCs.
%
% Inputs:
%    Properties possible to alter:
%
%    cellLocation            - location of bipolar RF center
%    patchSize               - size of retinal patch from sensor
%    timeStep                - time step of simulation from sensor
%    sRFcenter               - spatial RF of the center on receptor grid
%    sRFsurround             - spatial RF of the surround on receptor grid
%    temporalDifferentiator  - differentiator function
%    responseCenter          - Store the linear response of the center
%                              after convolution
%    responseSurround        - Store the linear response of the surround
%                              after convolution
%    bipolarResponseCenter   - store the bipolar response of the center
%                              after convolution
%    bipolarResponseSurround - store the bipolar response of the surround
%                              after convolution
%    tIR                     - bipolar temporal impulse response
%    threshold               - threshold placeholder data?
% 
% Outputs:
%    None
%

%% History 
%    05/xx/16  JRG (c) isetbio team
%    10/18/17  jnm  Comments & Formatting

% Examples:
%{
   cMosaic = coneMosaic; 
   bp = bipolarMosaic(cMosaic, 'on midget');
   bp.bipolarSet('sRFcenter',[0 0 0; 0 1 0; 0 0 0]);
   bp.bipolarSet('sRFsurround',[0 0 0; 0 1 0; 0 0 0]);
%}
%% Initialize and begin
narginchk(0, Inf);
p = inputParser; p.CaseSensitive = false; p.FunctionName = mfilename;
p.KeepUnmatched = true;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFieldsToSet = {...
    'cellLocation', ...
    'patchSize', ...
    'sRFcenter', ...
    'sRFsurround', ...
    'responseCenter', 'bipolarresponsecenter', ...
    'responseSurround', 'bipolarresponsesurround', ...
    'temporalDifferentiator', ...
    'threshold', ...
    'timeStep', ...
    'tIR'...
    };
p.addRequired('what', @(x) any(validatestring(ieParamFormat(x), ...
    allowableFieldsToSet)));
p.addRequired('value');
%% Parse and put results into structure p.
p.parse(varargin{:}); params = p.Results;

switch ieParamFormat(params.what)  % Lower case and remove spaces

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
        
    case{'threshold'}
        obj.threshold = params.value;
        
    case{'responsecenter'}
        obj.responseCenter = params.value;
        
    case{'responsesurround'}
        obj.responseSurround = params.value;
        
end
