function spatialArray(rgcM, varargin)
% Create RF centers and surrounds for an RGC mosaic
%
% Syntax:
%   spatialArray(rgcM, [varargin])
%
% Description:
%    The spatial RFs are generated according to the size of the pixel, cone
%    or bipolar mosaic, their spacing (in microns) and the diameter of the
%    RGC RF as determined by the TEE of the retinal patch.
%
%    The created values are then stored inside the provided rgc Mosaic,
%    with no individual values returned.
%
% Inputs:
%    rgcM                    - Object. A rgc Mosaic object.
%
% Outputs:
%   None.
%
% Optional key/value pairs:
%    centerSurroundSizeRatio - Numeric. The center surround size ratio.
%                              (Pass-through parameter from rgcRFEllipses)
%                              Default is sqrt(0.75)
%    centerSurroundAmpRatio  - Numeric. The center surround amplitude
%                              ratio. (Pass-through parameter from
%                              rgcRFEllipses) Default is 0.6 * 0.774.
%    centerNoise             - Numeric. The center noise. (Pass-through
%                              parameter from rgcRFEllipses) Default 0.15.
%    ellipseParams           - Matrix. A matrix of [A, B, rho]. (Pass-
%                              through parameter from rgcRFEllipses)
%                              Default is an empty matrix, [].
%    stride                  - Scalar. The stride length. Default [].
%
% Notes:
%    * See comments below for how the elliptical receptive fields are
%      created based on Chichilnisky and Kalmar logic (around line 60). The
%      code there and the ideas could be simplified over time.
%
% See Also:
%   initSpace (calls this function)
%

% History:
%    XX/XX/16  JRG/BW  (c) ISETBIO Team, 2016
%    06/12/19  JNM     Documentation pass

%% Manage parameters
p = inputParser;
p.KeepUnmatched = true;

p.addRequired('rgcM', @(x)(isa(rgcM, 'rgcMosaic')));
p.addParameter('stride', [], @isscalar);
p.parse(rgcM, varargin{:});

stride = p.Results.stride;
spread = rgcM.rfDiameter / 1;

% If stride is not set, use 2 std center to center separation.
if isempty(stride), stride = max(1, 1 * spread); end

%% Figure out the RGC cell locations with respect to the bipolar samples
[bpRow, bpCol, ~] = size(rgcM.input.cellLocation);

% rowSamples = 1:stride:bpRow;
% % The multiplicative factor on the col is hex spacing.
% colSamples = 1:((sqrt(3) / 2) * stride):bpCol;

rowSamples = 1:((sqrt(3) / 2) * stride):bpRow;
% The multiplicative factor on the col is hex spacing.
colSamples = 1:stride:bpCol;

rowSamples = rowSamples - mean(rowSamples);
colSamples = colSamples - mean(colSamples);

[Y, X] = meshgrid(rowSamples, colSamples);
cellLocation(:, :, 2) = Y;  % RGC cell center locations (Y)
cellLocation(:, :, 1) = X;  % RGC cell center locations (X)

%% Create elliptical spatial RFs for each cell
[rgcM.sRFcenter, rgcM.sRFsurround, rgcM.cellLocation, ...
    rgcM.ellipseMatrix] = rgcRFEllipses(cellLocation, spread, varargin{:});

end