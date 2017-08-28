function spatialArray(rgcM, varargin)
% BUILDSPATIALARRAY - Create RF centers and surrounds for an RGC mosaic
% 
% The spatial RFs are generated according to the size of the pixel, cone or
% bipolar mosaic, their spacing (in microns) and the diameter of the RGC RF
% as determined by the TEE of the retinal patch.
% 
%  [sRFcenter, sRFsurround, rfDiaMagnitude, ...
%    cellCenterLocations, tonicDrive, Qout] = 
%      buildSpatialRFArray(patchSizeMeters, nRowBipolars, nColBipolars, ...
%               rfDiameterMicrons, varargin)
% 
% Required inputs: 
%   patchSizeMeters   - Extent of the cone mosaic patch, inherited by
%                       bipolars and then here (meters)
%   nRowBipolars      - Number of row input samples from bpLayer
%   nColBipolars      - Number of col input samples from bpLayer
%   rfDiameterMicrons - Receptive field std in microns
%   
% Optional inputs:
%    centerNoise - Position jitter as a fraction of the bipolar size
%    baseLineFiringRate - Also called
%
% Outputs: 
%   sRFcenter           - intensity map as seen on the bipolar mosaic
%   sRFsurround         - intensity map as seen on the bipolar mosaic
%   cellCenterLocations - positions on the cone mosaic (microns)
%   tonicDrive          - tonic drive is the bias or DC term for the linear output of the GLM. 
%   Qout - ellipse parameters
%
% See also:  rgcInitSpace calls this function
% 
% Notes:  See comments below for how the elliptical receptive fields are
% created based on Chichilnisky and Kalmar logic (around line 60).  The
% code there and the ideas could be simplified over time.
%
% JRG/BW (c) ISETBIO Team, 2016

%% Manage parameters
p = inputParser;
p.KeepUnmatched = true;

p.addRequired('rgcM',@(x)(isa(rgcM,'rgcMosaic')));
p.addParameter('stride',[],@isscalar);
p.parse(rgcM,varargin{:});

stride = p.Results.stride;
spread = rgcM.rfDiameter/1;

% If stride is not set, use 2 std center to center separation.
if isempty(stride), stride = max(1,1*spread); end

%% Figure out the RGC cell locations with respect to the bipolar samples

%
[bpRow,bpCol,~] = size(rgcM.input.cellLocation);

% rowSamples = 1:stride:bpRow;   
% % The multiplicative factor on the col is hex spacing.
% colSamples = 1:((sqrt(3)/2)*stride):bpCol;

rowSamples = 1:((sqrt(3)/2)*stride):bpRow;   
% The multiplicative factor on the col is hex spacing.
colSamples = 1:stride:bpCol;

rowSamples = rowSamples - mean(rowSamples);
colSamples = colSamples - mean(colSamples);

[Y,X] = meshgrid(rowSamples,colSamples);
cellLocation(:,:,2) = Y; cellLocation(:,:,1) = X;   % RGC cell center locations

%% Create elliptical spatial RFs for each cell

[rgcM.sRFcenter, rgcM.sRFsurround, rgcM.cellLocation, rgcM.ellipseMatrix] = rgcRFEllipses(cellLocation, spread, varargin{:});

end