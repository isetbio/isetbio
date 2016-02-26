function obj = irSet(obj, varargin)
% rgcSet: a method of @rgc that sets rgc object 
% parameters using the input parser structure.
% 
%       val = rgcSet(rgc, property)
% 
% Inputs: rgc object, property to be set, value of property to be set
% 
% Outputs: object with property set
% 
% Proeprties:
%         name: type of rgc object, e.g., 'macaque RGC'
%         input: 'cone current' or 'scene RGB', depends on type of outer
%           segment object created.
%         temporalEquivEcc: the temporal equivalent eccentricity, used to 
%             determine the size of spatial receptive fields.   
%         mosaic: contains rgcMosaic objects for the five most common types
%           of RGCs: onParasol, offParasol, onMidget, offMidget,
%           smallBistratified.
%         numberTrials: the number of trials for spiking models LNP and GLM
% 
% Example:
%   rgc1 = rgcSet(rgc1, 'name', 'macaque RGC')
%   rgc1 = rgcSet(rgc1, 'temporalEquivEcc', 5)
% 
% 9/2015 JRG 

% Check for the number of arguments and create parser object.
% Parse key-value pairs.
% 

% % % We could do set using the superclass method
% obj = mosaicSet@rgcMosaic(obj, varargin{:});

% Check key names with a case-insensitive string, errors in this code are
% attributed to this function and not the parser object.
error(nargchk(0, Inf, nargin));
p = inputParser; p.CaseSensitive = false; p.FunctionName = mfilename;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFieldsToSet = {...   
        'name',...
        'input',...
        'temporalEquivEcc',...       
        'mosaic',...
        'numberTrials'...
    };
p.addRequired('what',@(x) any(validatestring(x,allowableFieldsToSet)));
p.addRequired('value');

% % Define what units are allowable.
% allowableUnitStrings = {'a', 'ma', 'ua', 'na', 'pa'}; % amps to picoamps
% 
% % Set up key value pairs.
% % Defaults units:
% p.addParameter('units','pa',@(x) any(validatestring(x,allowableUnitStrings)));

% Parse and put results into structure p.
p.parse(varargin{:}); params = p.Results;

% % Old error check on input.
% if ~exist('params','var') || isempty(params)
%     error('Parameter field required.');
% end
% if ~exist('val','var'),   error('Value field required.'); end;

% Set key-value pairs.
switch lower(params.what)
        
    case{'name'}
        obj.name = params.value;
    case{'input'}
        obj.input = params.value;
    case{'temporalequivecc'}        
        obj.temporalEquivEcc = params.value;
    case{'mosaic'}
        
        mosaicInd = length(obj.mosaic);
        if mosaicInd == 1 && isempty(obj.mosaic{mosaicInd})
            mosaicInd = 0;
        elseif mosaicInd >= 5
            mosaicInd = 0;
        end
        obj.mosaic{mosaicInd+1,1} = params.value;   
        
    case{'numbertrials'}
        cellTypes = length(obj.mosaic);
        if isa(obj.mosaic{1},'rgcLNP') | isa(obj.mosaic{1},'rgcGLM')| isa(obj.mosaic{1},'rgcPhys')
            
            for cellTypeInd = 1%:cellTypes
                obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd}, 'numberTrials', params.value);
            end
        else
            warning('The numberTrials property can only be set for rgcLNP, rgcGLM and rgcPhys models.');
        end
        
end

