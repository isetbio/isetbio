function obj = osSet(obj, varargin)
% Sets the isetbio outersegment object parameters.
%
% Syntax:
%   obj = osSet(obj, varargin)
%
% Description:
%    Set the isetbio outer segment object properties for the base class.
%
%    Examples are contained in the code. To access, type 'edit osSet.m'
%    into the Command Window.
%
% Inputs:
%    obj      - The isetbio outersegment object
%    varargin - The additional parameters required to assign a value to an
%               outersegment object parameter. Includes a parameter and a
%               value. Parameter options inclide:
%       {'noiseFlag'} - 'random'-Default, 'frozen', 'none'.
%       {'patchSize'} - Cone current as a function of time
%       {'timeStep'}  - Noisy cone current signal
%       {'size'}      - Array size of photon rate
%       {'coneCurrentSignal'}
%                     - Cone current as a function of time
%       {'sigma'}     - Rhodopsin activity decay rate (1/sec). Default 22
%       {'eta'}       - Phosphodiesterase activation rate constant (1/sec)
%                       Default 2000
%       {'beta'}      - Rate constant for calcium removal in 1/sec
%                       Default 9
%       {'opsingain'} - So stimulus can be in R*/sec (this is rate of
%                       increase in opsin activity per R*/sec). Default 10
%
% Outputs:
%    obj      - The modified isetbio outersegment object.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: JNM - Changed noiseFlag from 0/1 to none/frozen/random to
%      match with everywhere else.]
%

% Examples:
%{
    adaptedOS = osCreate;
	osGet(adaptedOS, 'noiseFlag')
    adaptedOS = osSet(adaptedOS, 'noiseFlag', 'frozen')
    adaptedOS = adaptedOS.osSet('noiseFlag', 'none')
%}

%%
% Check for the number of arguments and create parser object.
% Parse key-value pairs.
% 
% Check key names with a case-insensitive string, errors in this code are
% attributed to this function and not the parser object.
narginchk(0, Inf);
p = inputParser;
p.CaseSensitive = false;
p.FunctionName = mfilename;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFieldsToSet = {...
    'noiseflag', ...
    'conecurrentsignal', ...
    'patchsize', ...
    'timestep', ...
    'sigma', 'model.sigma', ...
    'eta', 'model.eta', ...
    'beta', 'model.beta', ...
    'opsingain', 'model.opsingain'};

p.addRequired('what', @(x) any(validatestring(ieParamFormat(x), ...
    allowableFieldsToSet)));
p.addRequired('value');

% Parse and put results into structure p.
p.parse(varargin{:}); params = p.Results;
%% Loop through param/value pairs

switch ieParamFormat(params.what)  % Lower case and remove space             
    case {'noiseflag'}
        obj.noiseFlag = params.value;

    case{'patchsize'}
        obj.patchSize = params.value;

    case{'timestep'}
        obj.timeStep = params.value;

    case{'conecurrentsignal'}
        obj.coneCurrentSignal = params.value;

    case{'sigma', 'model.sigma'}        
        % sigma = 10;
        % rhodopsin activity decay rate (1/sec) - default 22
        fprintf(['Choose a value within this range:\n sigma in ' ...
            'periphery is 22, sigma in fovea is 10']);
        obj.model.sigma = params.value;

    case{'eta', 'model.eta'}        
        % eta = 700;
        % phosphodiesterase activation rate constant (1/sec) - default 2000
        fprintf(['Choose a value within this range:\n eta in periphery' ...
            ' is 2000, eta in fovea is 700']);
        obj.model.eta = params.value;

    case{'beta', 'model.beta'}        
        % beta = 5;
        % rate constant for calcium removal in 1/sec - default 9
        fprintf(['Choose a value within this range:\n beta in ' ...
            'periphery is 9, beta in fovea is 5']);
        obj.model.beta = params.value;

    case{'opsingain', 'model.opsingain'}
        % OpsinGain = 12;
        % so stimulus can be in R*/sec (this is rate of increase in opsin
        % activity per R*/sec) - default 10
        fprintf(['Choose a value within this range:\n opsingain in ' ...
            'periphery is 10, opsingain in fovea is 12']);
        obj.model.OpsinGain = params.value;

end
