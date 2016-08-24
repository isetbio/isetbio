function val = osGet(obj, varargin)
% Gets the isetbio outersegment object parameters.
% 
% Parameters:
%       {'noiseFlag'} -  sets current as noise-free ('0') or noisy ('1')
%       {'patchSize'} - cone current as a function of time
%       {'timeStep'} - noisy cone current signal
%       {'size'} - array size of photon rate
%       {'coneCurrentSignal'} - cone current as a function of time
% 
% Model Parameters:
%       obj.model.sigma = 22;  % rhodopsin activity decay rate (1/sec) - default 22
%       obj.model.phi = 22;     % phosphodiesterase activity decay rate (1/sec) - default 22
%       obj.model.eta = 2000;	  % phosphodiesterase activation rate constant (1/sec) - default 2000
%       obj.model.gdark = 20.5; % concentration of cGMP in darkness - default 20.5
%       obj.model.k = 0.02;     % constant relating cGMP to current - default 0.02
%       obj.model.h = 3;       % cooperativity for cGMP->current - default 3
%       obj.model.cdark = 1;  % dark calcium concentration - default 1
%       obj.model.beta = 9;	  % rate constant for calcium removal in 1/sec - default 9
%       obj.model.betaSlow = 0.4; % rate constant for slow calcium modulation of channels - default 0.4
%       obj.model.n = 4;  	  % cooperativity for cyclase, hill coef - default 4
%       obj.model.kGc = 0.5;   % hill affinity for cyclase - default 0.5
%       obj.model.OpsinGain = 10; % so stimulus can be in R*/sec (this is rate of increase in opsin activity per R*/sec) - default 10
%        
%                     
% osGet(adaptedOS, 'noiseFlag')
% 
% 8/2015 JRG NC DHB


%% Check for the number of arguments and create parser object.
% 
% Check key names with a case-insensitive string, errors in this code are
% attributed to this function and not the parser object.
narginchk(0, Inf);
p = inputParser; p.CaseSensitive = false; p.FunctionName = mfilename;

% Make key properties that can be set required arguments.
allowableFieldsToSet = {...
    'noiseflag',...
    'patchsize',...
    'timestep',...
    'size',...
    'conecurrentsignal',...
    'model'...
    'sigma','model.sigma',...
    'eta','model.eta',...
    'beta','model.beta',...
    'opsingain','model.opsingain',...
    'gdark','model.gdark',...
    'k','model.k',...
    'h','model.h',...
    'cdark','model.cdark',...
    'betaslow','model.betaslow',...
    'n','model.n',...
    'kgc','model.kgc'};
p.addRequired('what',@(x) any(validatestring(ieParamFormat(x),allowableFieldsToSet)));

% Parse and put results into structure p.
p.parse(varargin{:}); params = p.Results;

switch ieParamFormat(params.what);  % Lower case and remove spaces

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
        
    case{'sigma','model.sigma'}
        % sigma = 10;  % rhodopsin activity decay rate (1/sec) - default 22
        val = obj.model.sigma;
        
    case{'eta','model.eta'}
        % eta = 700;      % phosphodiesterase activation rate constant (1/sec) - default 2000
        val = obj.model.eta;
        
    case{'beta','model.beta'}
        % beta = 5;      % rate constant for calcium removal in 1/sec - default 9
        val = obj.model.beta;
        
    case{'opsingain','model.opsingain'}
        % OpsinGain = 12; % so stimulus can be in R*/sec (this is rate of increase in opsin activity per R*/sec) - default 10
        val = obj.model.OpsinGain;
        
    case{'phi','model.phi'}
        val = obj.model.phi;
        
     case{'gdark','model.gdark'}
         val = obj.model.gdark;
         
    case{'k','model.k'}
        val = obj.model.k;
        
    case{'h','model.h'}
        val = obj.model.h;
        
    case{'cdark','model.cdark'}
        val = obj.model.cdark;
        
    case{'betaslow','model.betaslow'}
        val = obj.model.betaSlow;
        
    case{'n','model.n'}
        val = obj.model.n;
        
    case{'kgc','model.kgc'}
        val = obj.model.kGc;   
                                   
end

