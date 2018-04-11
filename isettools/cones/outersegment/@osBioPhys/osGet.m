function val = osGet(obj, varargin)
% Gets the isetbio outersegment object parameters.
%
% Syntax:
%    val = osGet(obj, varargin)
%
% Description:
%    Retrieve the isetbio outer segment object parameters.
%
%    Model Parameters:
%       All parameters are of the format obj.model.<param>, where param can
%       be any of the following:
%         sigma
%             - Rhodopsin activity decay rate (1/sec)
%             - Default 22
%         phi
%             - Phosphodiesterase activity decay rate (1/sec)
%             - Default 22
%         eta
%             - Phosphodiesterase activation rate constant (1/sec)
%             - Default 2000
%         gdark
%             - Concentration of cGMP in darkness
%             - Default 20.5
%         k
%             - Constant relating cGMP to current
%             - Default 0.02
%         h
%             - Cooperativity for cGMP->current
%             - Default 3
%         cdark
%             - Dark calcium concentration
%             - Default 1
%         beta
%             - Rate constant for calcium removal in 1/sec
%             - Default 9
%         betaSlow
%             - Rate constant for slow calcium modulation of channels
%             - Default 0.4
%         n
%             - Cooperativity for cyclase, hill coef
%             - Default 4
%         kGc
%             - Hill affinity for cyclase
%             - Default 0.5
%         OpsinGain
%             - So stimulus can be in R*/sec (this is rate of increase in
%               opsin activity per R*/sec)
%             - Default 10
%
%    There are examples in the fole. To access, type 'edit osGet.m' into
%    the Command Window.
%
% Inputs:
%    obj      - The isetbio outersegment object
%    varargin - The additional information, including the parameter to
%               retireve. Possible parameters include:
%         {'noiseFlag'} -  sets current as noise-free ('0') or noisy ('1')
%         {'patchSize'} - cone current as a function of time
%         {'timeStep'} - noisy cone current signal
%         {'size'} - array size of photon rate
%         {'coneCurrentSignal'} - cone current as a function of time
%
% Outputs:
%    val      - The value of the requested parameter.
%
% Optional key/value pairs:
%    None.
%

% History:
%    08/xx/15  JRG NC DHB  Created
%    02/15/18  jnm         Formatting
%    04/07/18  dhb         Fix example.

% Examples:
%{
    adaptedOS = osCreate;
    osGet(adaptedOS, 'noiseFlag')
%}

%% Check for the number of arguments and create parser object.
% 
% Check key names with a case-insensitive string, errors in this code are
% attributed to this function and not the parser object.
narginchk(0, Inf);
p = inputParser;
p.CaseSensitive = false;
p.FunctionName = mfilename;

% Make key properties that can be set required arguments.
allowableFieldsToSet = {...
    'noiseflag', ...
    'patchsize', ...
    'timestep', ...
    'size', ...
    'conecurrentsignal', ...
    'model'...
    'sigma', 'model.sigma', ...
    'eta', 'model.eta', ...
    'beta', 'model.beta', ...
    'opsingain', 'model.opsingain', ...
    'gdark', 'model.gdark', ...
    'k', 'model.k', ...
    'h', 'model.h', ...
    'cdark', 'model.cdark', ...
    'betaslow', 'model.betaslow', ...
    'n', 'model.n', ...
    'kgc', 'model.kgc'};
p.addRequired('what', @(x) any(validatestring(ieParamFormat(x), ...
    allowableFieldsToSet)));

% Parse and put results into structure p.
p.parse(varargin{:}); params = p.Results;

switch ieParamFormat(params.what)  % Lower case and remove spaces
    case {'noiseflag'}        
        val = obj.noiseFlag;

    case{'patchsize'}
        val = obj.patchSize;

    case{'timestep'}
        val = obj.timeStep;

    case{'size'}
        val = size(obj.coneCurrentSignal);

    case{'conecurrentsignal'}
        val = obj.coneCurrentSignal;

    case{'model'}
        val = obj.model;

    case{'sigma', 'model.sigma'}
        % sigma = 10;
        % rhodopsin activity decay rate (1/sec) - default 22
        val = obj.model.sigma;

    case{'eta', 'model.eta'}
        % eta = 700;
        % phosphodiesterase activation rate constant (1/sec) - default 2000
        val = obj.model.eta;

    case{'beta', 'model.beta'}
        % beta = 5;
        % rate constant for calcium removal in 1/sec - default 9
        val = obj.model.beta;

    case{'opsingain', 'model.opsingain'}
        % OpsinGain = 12;
        % so stimulus can be in R*/sec (this is rate of increase in opsin
        % activity per R*/sec) - default 10
        val = obj.model.OpsinGain;

    case{'phi', 'model.phi'}
        val = obj.model.phi;

     case{'gdark', 'model.gdark'}
         val = obj.model.gdark;

    case{'k', 'model.k'}
        val = obj.model.k;

    case{'h', 'model.h'}
        val = obj.model.h;

    case{'cdark', 'model.cdark'}
        val = obj.model.cdark;

    case{'betaslow', 'model.betaslow'}
        val = obj.model.betaSlow;

    case{'n', 'model.n'}
        val = obj.model.n;

    case{'kgc', 'model.kgc'}
        val = obj.model.kGc;   

end
