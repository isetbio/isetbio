function val = mosaicGet(obj, param, varargin)
% @rgcMosaic super class mosaicGet method
% 
%    val = mosaicGet(rgc.mosaic, property)
% 
% The subclasses have special values, and if the request is not for one of
% those then this call is invoked to get the value
%
% Inputs: 
%   obj    - rgc object
%   param  - parameter string
%   varargin - Not used yet, but will be used for units and other things.
% 
% Outputs: 
%   val - parameter value
% 
% Properties that can be gotten:
%    'cellType'         - type of RGC of which mosaic is composed
%    'rfDiameter'       - 1 stdev diameter in pixels of spatial RF
%    'rfDiaMagnitude'   - magnitude of spatial RF at 1 stdev
%    'cellLocation'     - coordinates of spatial RF center in ??? frame
%    'sRFcenter'        - center spatial RF surfaces
%    'sRFsurround'      - surround spatial RF surfaces
%    'tCenter'          - center temporal impulse response
%    'tSurround'        - surround temopral impulse response
%    'linearResponse'   - linear response of all cells
%    'mosaic size'      - Row/col for cells
%
% Examples:
%   val = mosaicGet(rgc1.mosaic{1}, 'cellType')
%   val = mosaicGet(rgc1.mosaic{3}, 'linearResponse')
% 
% 9/2015 JRG 

%% Parse
p = inputParser; 
p.CaseSensitive = false; 
p.FunctionName = mfilename;

% This flag causes the parser not to throw an error here in the superclass
% call. The subclass call will throw an error.
p.KeepUnmatched = true;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowFields = {...
        'celltype',...
        'rfdiameter',...
        'rfdiamagnitude',...
        'celllocation',...
        'srfcenter',...
        'srfsurround',...
        'tcenter',...
        'tsurround',...
        'linearresponse'...
        'mosaicsize'
    };
p.addRequired('param',@(x) any(validatestring(ieParamFormat(x),allowFields)));

% Parse and put results into structure p.
p.parse(param,varargin{:}); 
param = p.Results.param;

% @JRG - Please comment on the units
% Get 
switch ieParamFormat(param)
    
    case{'celltype'}
        val = obj.cellType;
    case{'rfdiameter'}
        val = obj.rfDiameter;
    case{'rfdiamagnitude'}
        val = obj.rfDiaMagnitude;
    case{'celllocation'}
        % @JRG:  Units are ?
        val = obj.cellLocation;
    case{'srfcenter'}
        val = obj.sRFcenter;
    case{'srfsurround'}
        val = obj.sRFsurround;
    case{'tcenter'}
        val = obj.tCenter;
    case{'tsurround'}
        val = obj.tSurround;
    case{'linearresponse'}
        val = obj.linearResponse;
    case {'mosaicsize'}
        val = size(obj.cellLocation);
end

end
