function val = mosaicGet(obj, varargin)
% rgcMosaicGet: a method of @rgcMosaic that gets rgcMosaic object 
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
    'spikeResponse'...

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
    
    case{'parent'}
        val = obj.parent;        
    case{'namecelltype'}        
        val = obj.nameCellType;
    case{'receptivefielddiameter1std'}
        val = obj.receptiveFieldDiameter1STD;
    case{'spatialrfarray'}
        val = obj.spatialRFArray;
    case{'spatialrfonedim'}
        val = obj.spatialRFonedim;   
    case{'spatialrfcenter'}
        val = obj.spatialRFcenter;
    case{'spatialrfsurround'}
        val = obj.spatialRFsurround;
    case{'spatialrfcontours'}
        val = obj.spatialRFcontours;
    case{'spatialrffill'}
        val = obj.spatialRFFill;
    case{'cellcenterlocations'}
        val = obj.cellCenterLocations;
    case{'temporalimpulseresponse'}
        val = obj.temporalImpulseResponse;
    case{'linearresponse'}
        val = obj.linearResponse;
    case{'generatorfunction'}
        val = obj.generatorFunction;        
    case{'nlresponse'}        
        val = obj.nlResponse;
    case{'spikeresponse'}
        val = obj.spikeResponse;
        
end

