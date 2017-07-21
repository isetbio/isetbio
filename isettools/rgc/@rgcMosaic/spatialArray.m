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

spread = rgcM.rfDiameter/2;
if isempty(stride), stride = max(1,round(spread)); end

[nRowBipolars, nColBipolars,~] = size(rgcM.input.cellLocation);

% Figure out how many RGCs we need
nRGC = zeros(2,1);
nRGC(1) = length(1:stride:nRowBipolars);
nRGC(2) = length(1:stride:nColBipolars);
nRGC(2) = floor((2/sqrt(3))*nRGC(2));   % JRG's magical hex

%% Converting the spatial units from microns to bipolar samples

% % p[atchSizeMeters is the width (columns) of the patch of cones sampled by
% % the bipolars. This arrives in meters, we convert row/col in microns
% % Row/Col in um
% patchRowColMicrons = patchSizeMeters*[(nRowBipolars/nColBipolars), 1]*1e6; 
% 
% % This is the diameter of the cone mosaic sampled by each bipolar in units
% % of microns.
% bipolarDiameterMicrons = (patchRowColMicrons(2) / nColBipolars);
% 
% % The spatial mosaic dimensions are in units of bipolar spacing
% 
% % The RGC rfDiameter in terms of bipolar array samples.  So, one step is
% % one step in the bipolar array. From here on out, the spatial coordinate
% % system is with respect to the bipolar array.
% rfDiameterBipolars = rfDiameterMicrons / bipolarDiameterMicrons;
% 
% % Determine the number of RGC samples in the hexagonal mosaic
% % nRGC: RGC cells, row col. 
% nRGC    = floor(patchRowColMicrons ./ rfDiameterMicrons); % number of RGCs in h, v direction
% 
% % Adjust the scale factor to account for the hexagonal packing of the RGC
% % mosaic
% nRGC(2) = floor((2/sqrt(3))*nRGC(2));

%% Create elliptical spatial RFs for each cell

[rgcM.sRFcenter, rgcM.sRFsurround, rgcM.cellLocation, rgcM.ellipseMatrix] = rgcRFEllipses(nRGC, rgcM.rfDiameter, varargin{:});


end