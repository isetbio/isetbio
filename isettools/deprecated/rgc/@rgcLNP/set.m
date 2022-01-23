function obj = set(obj, param, val, varargin)
% Sets a rgcLNP object parameter
%
% Syntax:
%   @rgcLNP.set(param, value, [varargin])
%
% Description:
%    Passes along parameters not found here to @rgcMosaic.set
%
%    Please note that all superclass parameter options are also supported,
%    although they are not listed below.
%
%    This function contains examples of usage inline. To access these, type
%    'edit @rgcLNP/set.m' into the Command Window.
%
% Inputs:
%    obj   - Object. A rgc object.
%    param - String. The parameter string. Options, which include the
%            corresponding type for val, are:
%       tonicdrive: Cell. A Cell matrix, containing the baseline rate of
%                   the likelihood of spiking, which creates a nonzero
%                   firing rate in response to no input.
%       generatorfunction: Matrix. Convert the conditional intensity to the
%                          likelihood of observing a spike within a given
%                          time bin, used also in GLM case.
%       postspikefilter: Numeric. The post spike filter.
%       responsevoltage: Matrix. The nonlinear voltage response after
%                        application of the generator function and the
%                        spike coupling responses is represented here in
%                        the form of a matrix.
%    val   - VARIES. The parameter value. See param above for more details.
%
% Outputs:
%    obj   - Object. The obj with property set appropriately.
%
% Optional key/value pairs:
%    Needs to be added.
%

% History:
%    09/XX/15  JRG  (c) isetbio team
%    07/XX/16  JRG  updated
%    06/19/19  JNM  Documentation pass


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
p.addRequired('param', @ischar);
p.addRequired('val');

% Parse and put results into structure p.
p.parse(param, val, varargin{:});
param = ieParamFormat(p.Results.param);
val = p.Results.val;

%% Set key-value pairs
switch param
    % LNP subclass parameters
    case{'tonicdrive'}
        % The baseline rate of the likelihood of spiking, which creates a
        % nonzero firing rate in response to no input
        for ci1 = 1:size(obj.tonicDrive, 1)
            for ci2 = 1:size(obj.tonicDrive, 2)
                obj.tonicDrive{ci1, ci2} = params.value;
            end
        end
    case{'generatorfunction'}
        % Converts the conditional intensity to the likelihood of observing
        % a spike within a given time bin, used also in GLM case.
        obj.generatorFunction = params.value;
    case{'postspikefilter'}
        obj.postSpikeFilter = params.value;
    case{'responsevoltage'}
        % The nonlinear voltage response after application of the generator
        % function and the spike coupling responses is represented here
        obj.responseVoltage = params.value;
    otherwise
        % Super class parameters
        set@rgcMosaic(obj, param, val, varargin{:});
end

