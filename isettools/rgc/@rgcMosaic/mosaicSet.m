function obj = mosaicSet(obj, param, val, varargin)
% mosaicSet for @rgcMosaic superclass
%
%  mosaic = mosaicGet(mosaic, property, value, varargin)
%  
% The subclasses rgcLinear, rgcLNP and others use this mosaicSet unless
% they have a special variable.
%
% Inputs: 
%   mosaic - 
%   param  -
%   value  - 
% 
% Outputs: 
%    obj with property set appropriately
% 
% Properties that can be set:
%   'cellType',...        - type of RGC of which mosaic is composed
%   'rfDiameter',...      - 1 stdev diameter in pixels of spatial RF
%   'rfDiaMagnitude',...  - magnitude of spatial RF at 1 stdev
%   'cellLocation',...    - row col of center of spatial RF
%   'sRFcenter',...       - center spatial RF surfaces
%   'sRFsurround',...     - surround spatial RF surfaces
%   'tCenter',...         - center temporal impulse response
%   'tSurround',...       - surround temopral impulse response
%   'linearResponse',...  - linear response of all cells
% 
% Examples:
%   rgc1.mosaic{1} = mosaicSet(rgc1.mosaic{1}, 'cellType', 'onParasol')
%   rgc1.mosaic{1} = mosaicSet(rgc1.mosaic{1}, 'linearResponse', linearResponse)
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
        'rfdiaMagnitude',...
        'celllocation',...
        'srfcenter',...
        'srfsurround',...
        'tcenter',...
        'tsurround',...
        'tcenterall',...
        'tsurroundall',...
        'responselinear'...
        'responsespikes',...
        'dt'
    };
p.addRequired('param',@(x) any(validatestring(ieParamFormat(x),allowFields)));
p.addRequired('val');

p.parse(param,val,varargin{:}); 
param = p.Results.param;
val   = p.Results.val;

%% Set key-value pairs.
switch lower(param)
    case{'celltype'}
        obj.cellType = val;
    case{'rfdiameter'}
        obj.rfDiameter = val;
    case{'rfdiamagnitude'}
        obj.rfDiaMagnitude = val;
    case{'celllocation'}
        obj.cellLocation = val;
    case{'srfcenter'}
        obj.sRFcenter = val;
    case{'srfsurround'}
        obj.sRFsurround = val;
    case{'tcenter'}
        obj.tCenter = val;                
    case{'tsurround'}
        obj.tSurround = val;
                
    case{'tcenterall'}
        nCells = size(obj.sRFcenter);
        tCenterNew = cell(nCells(1),nCells(2));
        for ii = 1:nCells(1)
            for jj = 1:nCells(2)                               
                tCenterNew{ii,jj} = val;
            end
        end
        obj.tCenter = tCenterNew;
        
    case{'tsurroundall'}
        nCells = size(obj.sRFsurround);
        tSurroundNew = cell(nCells(1),nCells(2));
        for ii = 1:nCells(1)
            for jj = 1:nCells(2)                               
                tSurroundNew{ii,jj} = val;
            end
        end
        obj.tSurround = tSurroundNew;
    case {'dt'}
        obj.dt = val;
        
    case{'responselinear'}
        % Need to deal with possibility of multiple trials!
        obj.responseLinear = val;
    case {'responsespikes'}
        obj.responseSpikes = val;
end

end
