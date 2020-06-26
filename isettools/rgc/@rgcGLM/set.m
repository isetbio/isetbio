function obj = set(obj, param, val, varargin)
%  Sets a property for an rgcGLM object.
%
% Syntax:
%   obj = obj.set(param, val, [varargin])
%   obj = set(obj, param, val, [varargin])
%
% Description:
%    Set a property for a rgcGLM object.
%
%    This function contains examples of usage inline. To access these, type
%    'edit @rgcGLM/set.m' into the Command Window.
%
% Inputs:
%    obj    - Object. A rgc object.
%    param  - String. A parameter string, rgcGLM-specific options and their
%             corresponding value (val) type are:
%       generatorfunction:
%       numbertrials: Numeric. The number of trials for which spikes are
%                     computed from the linear response.
%       responsevoltage: Numeric. The "membrane voltage" from the spike
%                        computation in units of conditional intensity.
%                        This signal contains the effects of the post spike
%                        filter and coupling filters for individual spikes
%                        in a given trial.
%       postspikefilter: Numeric. The post spike filter, in units of
%                        conditional intensity.
%       couplingfilter: Numeric. The coupling filter in units of
%                       conditional intensity.
%       couplingmatrix: Matrix. The weights on the coupling filter between
%                       cells in normailzed units between zero and one.
%    val    - VARIES. The parameter value, see param for more information.
%
% Outputs:
%    obj - Object. The modified rgcGLM object.
%
% Optional key/value pairs:
%    Needs to be added.
%

% History:
%    09/XX/15  JRG  (c) isetbio team
%    07/XX/16  JRG  updated
%    06/21/19  JNM  Documentation pass

% Examples:
%{
    % ETTBSkip - skipping broken examples
    %
    % This needs code to define rgc1 before it could
    % possibly work.
    rgc1.mosaic{1} = mosaicSet(rgc1.mosaic{1}, 'cellType', 'onParasol')
    rgc1.mosaic{1} = mosaicSet(rgc1.mosaic{1}, 'psthResponse', psth)
%}

%% Parse
p = inputParser;
p.CaseSensitive = false;
p.FunctionName = mfilename;

p.addRequired('param');
p.addRequired('val');

p.parse(param, val, varargin{:});
param = ieParamFormat(p.Results.param);
val = p.Results.val;

switch param
    % GLM subclass parameters
    case{'generatorfunction'}
        % An inline function, e.g. @exp
        obj.generatorFunction = val;
    case{'numbertrials'}
        % The number of trials for which spikes are computed from the
        % linear response
        obj.numberTrials = val;
    case{'responsevoltage'}
        % The "membrane voltage" from the spike computation in units of
        % conditional intensity. This signal contains the effects of the
        % post spike filter and coupling filters for individual spikes in a
        % given trial.
        obj.responseVoltage = val;
    case{'postspikefilter'}
        % The post spike filter, in units of conditional intensity
        obj.postSpikeFilter = val;
    case{'couplingfilter'}
        % The coupling filter in units of conditional intensity.
        obj.couplingFilter = val;
    case{'couplingmatrix'}
        % The weights on the coupling filter between cells in normailzed
        % units between zero and one.
        obj.couplingMatrix = val;
    otherwise
        % Superclass
        set@rgcMosaic(obj, param, val, varargin{:});
end

end
