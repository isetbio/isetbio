function obj = mosaicSet(obj, varargin)
% rgcMosaicSet: a method of @rgcMosaic that sets rgcMosaic object 
% parameters using the input parser structure.
% 
% Parameters:
%       {''} -  
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
    'parent',...
    'input',...
    'nameCellType',...
    'receptiveFieldDiameter1STD',...
    'spatialRFArray',...
    'spatialRFonedim',...
    'spatialRFcenter',...
    'spatialRFsurround',...
    'spatialRFcontours',...
    'spatialRFFill',...
    'cellCenterLocations',...
    'temporalImpulseResponse',...
    'linearResponse',...
    'generatorFunction',...
    'nlResponse;',...
    'spikeResponse',...
    'couplingFilter',...
    'couplingMatrix',...
    'rasterResponse',...
    'psthResponse'...

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
    
    case{'parent'}
        obj.parent = params.value;           
    case{'input'}
        obj.input = params.value;
    case{'namecelltype'}        
        obj.nameCellType = params.value;
    case{'receptivefielddiameter1std'}
        obj.receptiveFieldDiameter1STD = params.value;
    case{'spatialrfarray'}
        obj.spatialRFArray = params.value;
    case{'spatialrfonedim'}
        obj.spatialRFonedim = params.value;
    case{'spatialrfcenter'}
        obj.spatialRFcenter = params.value;
    case{'spatialrfsurround'}
        obj.spatialRFsurround = params.value;
    case{'spatialrfcontours'}
        obj.spatialRFcontours = params.value;
    case{'spatialrffill'}
        obj.spatialRFFill = params.value;
    case{'cellcenterlocations'}
        obj.cellCenterLocations = params.value;
    case{'temporalimpulseresponse'}
        obj.temporalImpulseResponse = params.value;
    case{'linearresponse'}
        obj.linearResponse = params.value;
    case{'generatorfunction'}
        obj.generatorFunction = params.value;        
    case{'nlresponse'}        
        obj.nlResponse = params.value;
    case{'spikeresponse'}
        obj.spikeResponse = params.value;        
    case{'couplingfilter'}
        obj.couplingFilter = params.value;
    case{'couplingmatrix'}
        obj.couplingMatrix = params.value;
    case{'rasterresponse'}
        obj.rasterResponse = params.value;
    case{'psthresponse'}
        obj.psthResponse = params.value;
end

