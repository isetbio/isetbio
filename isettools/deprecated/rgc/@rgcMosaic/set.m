function obj = set(obj, param, val, varargin)
% Sets a property for an rgcMosaic object.
%
% Syntax:
%   mosaic = @rgcMosaic.set(mosaic, property, value, [varargin])
%
% Description:
%    The mosaicSet function sets a property for the mosaic object if the
%    property is defined in the rgcMosaic superclass. The subclasses of
%    rgcMosaic have properties not included here, and if the request is not
%    for one of the properties common to all subclasses, then the subclass
%    mosaicSet is called.
%
%    This function contains examples of usage inline. To access these, type
%    'edit @rgcMosaic\set' into the Command Window.
%
% Inputs:
%   obj    - Object. A rgc object.
%   param  - String. The parameter string. Some possible options are:
%         celltype: String. The RGC type
%         rfdiameter: Numeric. The spatial RF diameter in micrometers.
%         rfdiaMagnitude: Numeric. The magnitude of the linear spatial RF
%                         in units of spikes/sec at which one standard
%                         deviation contours are computed.
%         cellLocation: <Type>. The location of the RF center in units of
%                       zero-centered cone position, where the center of
%                       the RGC mosaic is at [0 0].
%         srfcenter: <Type>. Center of spatial RF surfaces.
%         srfsurround: <Type>. Surround of spatial RF surfaces.
%         tcenter: <Type>. The center temporal impulse response.
%         tsurround: <Type>. The surround temporal impulse response.
%         tcenterall: <Type>. <Fill in>.
%         tsurroundall: <Type>. <Fill in>.
%         responsevoltage: <Type>. Nonlinear voltage for GLM model
%         responselinear: <Type>. The linear response of all cells.
%         responsespikes: Array. The spike times in a given trial.
%         dt: Numeric. The bin subsampling size.
%   val    - VARIES. The parameter value, see param options for more
%            information about type and expected value(s).
%
% Outputs:
%    obj   - Object. The object with property set appropriately.
%
% Optional key/value pairs:
%    Needs to be added.
%
% Notes:
%    * See the settable responses on line 51.
%

% History:
%    09/XX/15  JRG  (c) isetbio team
%    07/XX/16  JRG  updated
%    06/12/19  JNM  Documentation pass

% Examples:
%{
    % ETTBSkip - skipping broken examples
    %
    % This needs code to define rgc1 before it could
    % possibly work.
    rgc1.mosaic{1} = mosaicSet(rgc1.mosaic{1}, 'cellType', 'onParasol')
    rgc1.mosaic{1} = mosaicSet(rgc1.mosaic{1}, ...
        'linearResponse', linearResponse)
%}

%% Parse
p = inputParser;
p.CaseSensitive = false;
p.FunctionName = mfilename;

% This flag causes the parser not to throw an error here in the superclass
% call. The subclass call will throw an error.
p.KeepUnmatched = true;

% [Note: BW - What are the spatial units?  On the mosaic of the bipolar
% input grid, I think.]
allowFields = {'celltype', 'rfdiameter', 'rfdiaMagnitude', ...
    'celllocation', 'srfcenter', 'srfsurround', 'tcenter', 'tsurround', ...
    'tcenterall', 'tsurroundall', 'responsevoltage', 'responselinear'...
    'responsespikes', 'dt'};
p.addRequired('param', ...
    @(x) any(validatestring(ieParamFormat(x), allowFields)));
p.addRequired('val');

p.parse(param, val, varargin{:});
param = ieParamFormat(p.Results.param);
val = p.Results.val;

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
        % If ir.mosaic{1}.cellLocation{1, 1} = [-40 -40], then the mosaic
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
        tCenterNew = cell(nCells(1), nCells(2));
        for ii = 1:nCells(1)
            for jj = 1:nCells(2)
                tCenterNew{ii, jj} = val;
            end
        end
        obj.tCenter = tCenterNew;
    case{'tsurroundall'}
        nCells = size(obj.sRFsurround);
        tSurroundNew = cell(nCells(1), nCells(2));
        for ii = 1:nCells(1)
            for jj = 1:nCells(2)
                tSurroundNew{ii, jj} = val;
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
