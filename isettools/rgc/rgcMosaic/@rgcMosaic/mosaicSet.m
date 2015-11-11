function obj = mosaicSet(obj, varargin)
% rgcMosaicSet: a method of @rgcMosaic that sets rgcMosaic object 
% parameters using the input parser structure.
% 
%       rgc.mosaic = mosaicGet(rgc.mosaic, property, value)
%  
% Inputs: rgc object, key-value pair of property and value to which it is
%   being set.
% 
% Outputs: obj with property set appropriately.% 
% 
% Properties that can be set:
%         'cellType',...        - type of RGC of which mosaic is composed
%         'rfDiameter',...      - 1 stdev diameter in pixels of spatial RF
%         'rfDiaMagnitude',...  - magnitude of spatial RF at 1 stdev
%         'cellLocation',...    - coordinates of center of spatial RF
%         'sRFcenter',...       - center spatial RF surfaces
%         'sRFsurround',...     - surround spatial RF surfaces
%         'tCenter',...         - center temporal impulse response
%         'tSurround',...       - surround temopral impulse response
%         'linearResponse',...  - linear response of all cells
% 
% 
% Examples:
%   rgc1.mosaic{1} = mosaicSet(rgc1.mosaic{1}, 'cellType', 'onParasol')
%   rgc1.mosaic{1} = mosaicSet(rgc1.mosaic{1}, 'linearResponse', linearResponse)
% 
% 9/2015 JRG  

% Check for the number of arguments and create parser object.
% Parse key-value pairs.
% 
% Check key names with a case-insensitive string, errors in this code are
% attributed to this function and not the parser object.
error(nargchk(0, Inf, nargin));
p = inputParser; p.CaseSensitive = false; p.FunctionName = mfilename;

% This flag causes the parser not to throw an error here in the superclass
% call. The subclass call will throw an error.
p.KeepUnmatched = true;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFieldsToSet = {...
        'cellType',...
        'rfDiameter',...
        'rfDiaMagnitude',...
        'cellLocation',...
        'sRFcenter',...
        'sRFsurround',...
        'tCenter',...
        'tSurround',...
        'linearResponse'...
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
    case{'celltype'}
        obj.cellType = params.value;
    case{'rfdiameter'}
        obj.rfDiameter = params.value;
    case{'rfdiamagnitude'}
        obj.rfDiaMagnitude = params.value;
    case{'celllocation'}
        obj.cellLocation = params.value;
    case{'srfcenter'}
        obj.sRFcenter = params.value;
    case{'srfsurround'}
        obj.sRFsurround = params.value;
    case{'tcenter'}
        obj.tCenter = params.value;
    case{'tsurround'}
        obj.tSurround = params.value;
    case{'linearresponse'}
        obj.linearResponse = params.value;
end

