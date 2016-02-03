function val = rgcGet(obj, varargin)
% rgcGet: a method of @rgc that gets rgc object 
% parameters using the input parser structure.
% 
%       val = rgcGet(rgc, property)
% 
% Inputs: rgc object, property to be gotten
% 
% Outputs: val of property
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
% 
% Example:
%   val = rgcGet(rgc1, 'name')
%   val = rgcGet(rgc1, 'input')
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
        'featureVector'...
%         'linearResponse'...
    };
p.addRequired('what',@(x) any(validatestring(x,allowableFieldsToSet)));

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
        val = obj.name;
    case{'input'}
        val = obj.input;
    case{'temporalequivecc'}        
        val = obj.temporalEquivEcc;
    case{'mosaic'}        
        val = obj.mosaic;                
    case{'featurevector'}
        val = [];
        cellTypes = length(obj.mosaic);
        for cellTypeInd = 1:cellTypes
            numberSpikes = mosaicGet(obj.mosaic{cellTypeInd}, 'rasterResponse');
            numberTrials = mosaicGet(obj.mosaic{cellTypeInd}, 'numberTrials');
            % obj.mosaic{1}.rasterResponse{1,1} is an array of NxTrials,
            % but since some of the array entries are zero, we count the
            % number of spikes per trial by finding the number of nonzero
            % array entries and take the mean.
            
            % Do this operation for each cell
            valC = (cellfun(@(x) mean(size(x, 1) - sum(isinf(1./x))),numberSpikes, 'un',0));
            
            % put into array
            valA = cell2mat(valC);
            val = [val; valA(:)];
            
        end
end

