function obj = set(obj, param, val, varargin)
% Sets a property for an rgcMosaic object.
%
%  mosaic = @rgcMosaic.set(mosaic, property, value, varargin)
%  
% The mosaicSet function sets a property for the mosaic object if the
% property is defined in the rgcMosaic superclass. The subclasses of
% rgcMosaic have properties not included here, and if the request is not
% for one of the properties common to all subclasses, then the subclass
% mosaicSet is called.
%
% Inputs: 
% 
%   obj    - rgc object
%   param  - parameter string
%   val    - parameter value
%   varargin - Not used yet, but will be used for units and other things.
% 
% Outputs: 
%    obj with property set appropriately
% 
% See the settable responses on line 51.
%
% Examples:
%   rgc1.mosaic{1} = mosaicSet(rgc1.mosaic{1}, 'cellType', 'onParasol')
%   rgc1.mosaic{1} = mosaicSet(rgc1.mosaic{1}, 'linearResponse', linearResponse)
% 
% 9/2015 JRG (c) isetbio team
% 7/2016 JRG updated

%% Parse
p = inputParser;
p.CaseSensitive = false; 
p.FunctionName = mfilename;

% This flag causes the parser not to throw an error here in the superclass
% call. The subclass call will throw an error.
p.KeepUnmatched = true;

% Make key properties that can be set required arguments, and require
% values along with key names.

%   'rfDiameter',...      - 1 stdev diameter in pixels of spatial RF
%   'rfDiaMagnitude',...  - magnitude of spatial RF at 1 stdev
%   'cellLocation',...    - row col of center of spatial RF
%   'sRFcenter',...       - center spatial RF surfaces
%   'sRFsurround',...     - surround spatial RF surfaces
%   'tCenter',...         - center temporal impulse response
%   'tSurround',...       - surround temopral impulse response
%   'linearResponse',...  - linear response of all cells
% 

% What are the spatial units?  On the mosaic of the bipolar input grid, I
% think (BW).
allowFields = {...
        'celltype',...        % RGC type
        'rfdiameter',...      % In samples of the bipolar mosaic ?
        'rfdiaMagnitude',...  % Not sure
        'celllocation',...    % In samples of the bipolar mosaic?
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
param = ieParamFormat(p.Results.param);
val   = p.Results.val;

%% Set key-value pairs.
switch lower(param)
    case{'celltype'}
        % String that stores cell type name
        obj.cellType = val;
    case{'rfdiameter'}
        % Spatial RF diameter in micrometers
        obj.rfDiameter = val;
    case{'rfdiamagnitude'}
        % Magnitude of linear spatial RF in units of spikes/sec at which 1
        % SD contours are computed.
        obj.rfDiaMagnitude = val;
    case{'celllocation'}        
        % Location of RF center in units of zero-centered cone position
        % The center of the RGC mosaic is at [0 0]
        % If ir.mosaic{1}.cellLocation{1,1} = [-40 -40], then the mosaic
        % underlying cone mosaic is about [80x80] depending on RF size 
        obj.cellLocation = val;
    case{'srfcenter'}
        % Linear spatial center RF in units of conditional intensity,
        % related by Poisson firing to spikes/sec.
        obj.sRFcenter = val;
    case{'srfsurround'}
        % Linear spatial surround RF in units of conditional intensity,
        % related by Poisson firing to spikes/sec.spikes/sec.
        obj.sRFsurround = val;
    case{'tcenter'}
        % Linear temporal center impulse response in units of conditional
        % intensity, related by Poisson firing to spikes/sec
        obj.tCenter = val;                
    case{'tsurround'}
        % Linear temporal surround impulse response in units of conditional
        % intensity, related by Poisson firing to spikes/sec
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
        % The bin subsampling size. In the original Pillow code, was a
        % fraction of the sampling rate of the linear response (usually
        % 1/120 = 0.0083 sec). Now it takes into account the linear
        % sampling rate and is given in units of microseconds.
        obj.dt = val;
        
    case{'responselinear'}
        % Linear response in units of conditional intensity, related by
        % Poisson firing to spikes/sec
        obj.responseLinear = val;
    case {'responsespikes'}
        % The spike times on a given trial.
        obj.responseSpikes = val;
end

end
