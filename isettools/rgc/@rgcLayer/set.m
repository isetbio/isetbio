function obj = set(obj, param, val, varargin)
% A method of @rgc that sets rgc object parameters using an input parser.
%
% syntax:
%   val = @rgcLayer.set(ir, param, value, varargin)
%
% Description:
%    This function is a method of @rgc that sets rgc object parameters
%    using an input parser structure.
%
%    This function contains examples of usage inline. To access these, type
%    'edit @rgcLayer/set.m' into the Command Window.
%
% Inputs:
%    obj     - Object. A RGC layer object
%    param   - String. The parameter whose value you wish to assign.
%              Options include the following:
%       name: String. A type of rgc object, e.g., 'macaque RGC'
%       input: String. The input type, such as 'cone current' or 'scene
%              RGB', which depends on type of outer segment object created.
%       temporalEquivEcc: Numeric. The temporal equivalent eccentricity,
%                         which is used to determine the size of spatial
%                         receptive fields.
%       timing: Numeric. The stimulus input time step.
%       spacing: Numeric. The space between mosaics.
%       mosaic: Object. This contains rgcMosaic objects for the five most
%               common types of RGCs: onParasol, offParasol, onMidget,
%               offMidget, and smallBistratified.
%       numberTrials: Numeric. The number of trials for spiking models LNP
%                     and GLM.
%       recordedSpikes: Cell. A physical record of the spikes for ISETBio.
%       numbertrials/ntrials: Numeric. The number of trials. [Note: XXX -
%                             At some point in the future this will be a
%                             vector containing the number of trials for
%                             each mosaic in the inner retina.]
%    val     - VARIES. The parameter to assign to the
%
% Outputs:
%    obj     - Object. The modified rgc layer object.
%
% Optional key/value pairs:
%    nMosaic - Numeric. The mosaic number. Default [].
%

% History:
%    09/XX/15  JRG  Created
%    06/20/19  JNM  Documentation pass

% Examples:
%{
    % ETTBSkip - skipping broken examples
    %
    % This needs code to define rgc1 before it could
    % possibly work.
    rgc1 = rgcSet(rgc1, 'name', 'macaque RGC')
    rgc1 = rgcSet(rgc1, 'temporalEquivEcc', 5)
%}

%% Parse
p = inputParser;
p.CaseSensitive = false;
p.FunctionName = mfilename;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFieldsToSet = {'name', 'input', 'temporalEquivEcc', 'timing', ...
    'spacing', 'mosaic', 'numberTrials', 'recordedSpikes', ...
    'numbertrials', 'ntrials'};
p.addRequired('param', @(x) any(validatestring(x, allowableFieldsToSet)));
p.addRequired('val');
p.addParameter('nMosaic', [], @isscalar);

% Parse and put results into structure p.
p.parse(param, val, varargin{:});
param = ieParamFormat(p.Results.param);
val = p.Results.val;

%% Set key-value pairs.
switch param
    case{'name'}
        obj.name = val;
    case{'input'}
        obj.input = val;
    case{'temporalequivecc'}
        obj.temporalEquivEcc = val;
    case{'mosaic'}
        % set('mosaic', val, 'nMosaic', integer);
        nMosaic = p.Results.nMosaic;
        if isempty(nMosaic)
            obj.mosaic{end + 1} = val;
        else
            obj.mosaic{nMosaic} = val;
        end
    case{'numbertrials', 'ntrials'}
        % Will be a vector some day for number of trials for each mosaic in
        % the ir. Now, it is just a plain old number.
        obj.nTrials = val;
    case{'timing'}
        obj.timing = val;
    case{'spacing'}
        obj.spacing = val;
    case{'recordedspikes'}
        % For Phys data into ISETBIO
        % @JRG - Comment. Preallocate space.
        for cellind = 1:length(val)
            for iTrial = 1:length(val{1}.rasters.recorded)
                recorded_spiketimes{cellind, 1, iTrial} = ...
                    (val{cellind}.rasters.recorded{iTrial});
            end
        end

        if isa(obj.mosaic{1}, 'rgcPhys')
            obj.mosaic{1} = mosaicSet(obj.mosaic{1}, ...
                'responseSpikes', recorded_spiketimes);
        end
end

end
