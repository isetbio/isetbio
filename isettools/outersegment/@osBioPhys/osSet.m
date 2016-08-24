function obj = osSet(obj, varargin)
% Sets the isetbio outersegment object parameters.
% 
% Parameters:
%       {'noiseFlag'} -  sets current as noise-free ('0') or noisy ('1')
%       {'patchSize'} - cone current as a function of time
%       {'timeStep'} - noisy cone current signal
%       {'size'} - array size of photon rate
%       {'coneCurrentSignal'} - cone current as a function of time
%       {'sigma'} - rhodopsin activity decay rate (1/sec) - default 22
%       {'eta'} - phosphodiesterase activation rate constant (1/sec) - default 2000
%       {'beta'} - rate constant for calcium removal in 1/sec - default 9
%       {'opsingain'} - so stimulus can be in R*/sec (this is rate of increase in opsin activity per R*/sec) - default 10
% 
% osGet(adaptedOS, 'noiseFlag')

%%
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
    'conecurrentsignal',...
    'patchsize',...
    'timestep',...
    'sigma','model.sigma',...
    'eta','model.eta',...
    'beta','model.beta',...
    'opsingain','model.opsingain'};

p.addRequired('what',@(x) any(validatestring(ieParamFormat(x),allowableFieldsToSet)));
p.addRequired('value');

% Parse and put results into structure p.
p.parse(varargin{:}); params = p.Results;
%% Loop through param/value pairs

switch ieParamFormat(params.what);  % Lower case and remove space
                    
    case {'noiseflag'}
        obj.noiseFlag = params.value;
        
    case{'patchsize'}
        obj.patchSize = params.value;
        
    case{'timestep'}
        obj.timeStep = params.value;
        
    case{'conecurrentsignal'}
        obj.coneCurrentSignal = params.value;
        
    case{'sigma','model.sigma'}        
        % sigma = 10;  % rhodopsin activity decay rate (1/sec) - default 22
        fprintf('Choose a value within this range:\n sigma in periphery is 22, sigma in fovea is 10');
        obj.model.sigma = params.value;
        
    case{'eta','model.eta'}        
        % eta = 700;      % phosphodiesterase activation rate constant (1/sec) - default 2000
        fprintf('Choose a value within this range:\n eta in periphery is 2000, eta in fovea is 700');
        obj.model.eta = params.value;
        
    case{'beta','model.beta'}        
        % beta = 5;      % rate constant for calcium removal in 1/sec - default 9
        fprintf('Choose a value within this range:\n beta in periphery is 9, beta in fovea is 5');
        obj.model.beta = params.value;
        
    case{'opsingain','model.opsingain'}
        % OpsinGain = 12; % so stimulus can be in R*/sec (this is rate of increase in opsin activity per R*/sec) - default 10
        fprintf('Choose a value within this range:\n opsingain in periphery is 10, opsingain in fovea is 12');
        obj.model.OpsinGain = params.value;
        
end

