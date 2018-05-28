function obj = set(obj, varargin)
% Syntax:
%
%   bp.set('PARAM1',var1, 'PARAM2', var2, ...);
%
% Description:
%    Gets isetbio bipolar object parameters.
%
%    The bipolar object allows the simulated cone responses to be passed on
%    to the inner retina object and to approxiately maintain its impulse
%    response. This will allow us to run the nonlinear biophysical cone
%    outer segment model and pass its results on to the bipolar stage and
%    then RGCs.
%
% Inputs:
%    Properties amenable to being set.
%
%    cellLocation            - location of bipolar RF center
%    patchSize               - size of retinal patch from sensor
%    timeStep                - time step of simulation from sensor
%    sRFcenter               - spatial RF of the center on receptor grid
%    sRFsurround             - spatial RF of the surround on receptor grid
%    temporalDifferentiator  - differentiator function
%    responseCenter          - Store linear response of the center after
%                              convolution
%    responseSurround        - Store linear response of the surround after
%                              convolution
%    bipolarResponseCenter   - store the bipolar response of the center
%                              after convolution
%    bipolarResponseSurround - store the bipolar response of the surround
%                              after convultion
%    threshold               - threshold placeholder data?
%    tIR                     - bipolar temporal impulse response
%

%% History
% 5/2016 JRG (c) isetbio team
%
%    10/18/17  jnm  Comments & formatting

%% Examples:
%{
   % ETTBSkip - Example is broken. Variables undefined. Remove this line when fixed.
   innerRetinaSU.mosaic{1}.set('dt',1);
%}

%% Initialize and begin
narginchk(0, Inf);
p = inputParser;
p.CaseSensitive = false;
p.FunctionName = mfilename;
p.KeepUnmatched = true;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFieldsToSet = {...
    'cellLocation', ...
    'patchSize', ...
    'responseCenter', 'bipolarresponsecenter', ...
    'responseSurround', 'bipolarresponsesurround'...
    'sRFcenter', ...
    'sRFsurround', ...
    'temporalDifferentiator', ...
    'threshold', ...
    'timeStep', ...
    'tIR', ...
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

