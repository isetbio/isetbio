function val = get(obj, param, varargin)
% A method of @rgc that gets rgc object parameters using an input parser.
%
% Syntax:
%   val = @rgcLayer.get(rgc, param, [varargin])
%
% Description:
%    A @rgc method that retrieves rgc object parameters using an input
%    parser structure.
%
%    This function contains examples of usage inline. To access these, type
%    'edit @rgcLayer/get.m' into the Command Window.
%
% Inputs:
%    rgc   - Object. An rgc object.
%    param - String. A string indicating which parameter to retrieve.
%            Specified parameter affects the output type. Parameter options
%            include the following:
%        name: String. The rgc object type. e.g. 'macaque RGC'.
%        input: String. The outer segment object type, either 'cone
%               concurrent' or 'scene RGB'.
%        temporalEquivEcc: Numeric. The temporal Equivalent Eccentricity.
%                          The eccentricity (in mm) of the location in the
%                          nasal half is mapped onto the equivalent
%                          temporal eccentricity. This is equivalent with
%                          respect to the RGC density. Used to determine
%                          the size of the spatial receptive fields.
%        mosaic: Object. The rgc mosaic objects for the five most common
%                types of RGCs: onParasol, offParasol, onMidget, offMidget,
%                and smallBistratified.
%        timing: Numeric. The stimulus input time step
%        spacing: Numeric. The space between mosaics.
%        numbertrials: Numeric. The number of trials.
%
% Outputs:
%   val   - VARIES. The value of the property. See param above for possible
%           return types.
%
% Optional key/value pairs:
%    Needs to be added.
%

% History:
%    XX/XX/15  JRG/BW  ISETBIO Team, 2015
%    06/21/19  JNM     Documentation pass

% Examples:
%{
    % ETTBSkip - skipping broken examples
    %
    % This needs code to define rgc1 before it could
    % possibly work.
    val = rgcGet(rgc1, 'name')
    val = rgcGet(rgc1, 'input')
%}

%% Parse
p = inputParser;
p.CaseSensitive = false;
p.FunctionName = mfilename;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFieldsToSet = {'name', 'input', 'temporalEquivEcc', 'mosaic', ...
    'timing', 'spacing', 'numbertrials'};
p.addRequired('param', ...
    @(x) any(validatestring(ieParamFormat(x), allowableFieldsToSet)));

% Parse and put results into structure p.
p.parse(param, varargin{:});
param = ieParamFormat(p.Results.param);

%% Set key-value pairs.
switch param
    case{'name'}
        val = obj.name;
    case{'input'}
        val = obj.input;
    case{'temporalequivecc'}
        val = obj.temporalEquivEcc;
    case{'mosaic'}
        val = obj.mosaic;
    case{'timing'}
        val = obj.timing;
    case{'spacing'}
        val = obj.spacing;
    case {'numbertrials'}
        val = obj.nTrials;
end

end
