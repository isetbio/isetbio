function val = get(obj, param, varargin)
% Gets a property from an rgcGLM object.
%
% Syntax:
%   val = @rgcGLM.get(param, [varargin])
%
% Description:
%    Retrieve a rgcGLM property.
%
%    This function contains examples of usage inline. To acces these, type
%    'edit @rgcGLM/get' into the Command Window.
%
% Inputs:
%    obj    - Object. A rgc object.
%    param  - String. A parameter string. Determines the type of the output
%             val, with options such as:
%       generatorfunction: Matrix. Convert the conditional intensity to the
%                          likelihood of observing a spike within a given
%                          time bin, used also in GLM case.
%       postSpikeFilter: Numeric. The post spike filter, in units of
%                        conditional intensity.
%       numberTrials: Numeric. The number of trials for which spikes are
%                              computed from the linear response.
%       responseVoltage: Matrix. The nonlinear voltage response after
%                        application of the generator function and the
%                        spike coupling responses is represented here in
%                        the form of a matrix.
%       couplingFilter: Numeric. The coupling filter in units of
%                       conditional intensity.
%       couplingMatrix: Matrix. The weights on the coupling filter between
%                       cells in normailzed units between zero and one.
%
% Outputs:
%    val    - VARIES. The parameter value, see param above for details.
%
% Optional key/value pairs:
%    None.
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
    val = mosaicGet(rgc1.mosaic{1}, 'cell type')
    val = mosaicGet(rgc1.mosaic{3}, 'psth response')
%}

%% Parse
p = inputParser;
p.CaseSensitive = false;
p.FunctionName = mfilename;
p.KeepUnmatched = true;

p.addRequired('param');
p.parse(param, varargin{:});
param = ieParamFormat(p.Results.param);

%% Set key-value pairs.
switch param
    % Specific to the GLM case
    case{'generatorfunction'}
        % An inline function, e.g. @exp
        val = obj.generatorFunction;
    case{'postspikefilter'}
        % The post spike filter, in units of conditional intensity
        val = obj.postSpikeFilter;
    case{'numbertrials'}
        % The number of trials for which spikes are computed from the
        % linear response
        if ~isempty(obj.responseSpikes)
            val = size(obj.responseSpikes, 3);
        else
            val = 0;
        end
    case{'responsevoltage'}
        % The "membrane voltage" from the spike computation in units of
        % conditional intensity. This signal contains the effects of the
        % post spike filter and coupling filters for individual spikes in a
        % given trial.
        val = obj.responseVoltage;
    case{'couplingfilter'}
        % The coupling filter in units of conditional intensity.
        val = obj.couplingFilter;
    case{'couplingmatrix'}
        % The weights on the coupling filter between cells in normailzed
        % units between zero and one.
        val = obj.couplingMatrix;
    otherwise
        % Super class
        val = get@rgcMosaic(obj, param, varargin{:});
end

end
