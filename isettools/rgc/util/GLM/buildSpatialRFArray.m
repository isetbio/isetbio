function [sRFcenter, sRFsurround, cellCenterLocations, tonicDrive, Qout] = ...
    buildSpatialRFArray(patchSizeMeters, nRowBipolars, nColBipolars, rfDiameterMicrons, varargin)
% Builds the spatial RF center and surround for the cells in a mosaic
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
%   patchSizeMeters   - Center to center of the RF in microns
%   nRowBipolars      - Number of input samples
%   nColBipolars      - Number of input samples
%   rfDiameterMicrons - Receptive field of 1 std in microns
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
% Example:
%   Build spatial RFs of the RGCs in this mosaic
% 
% JRG (c) isetbio Team

%% Manage parameters
p = inputParser;
p.KeepUnmatched = true;

p.addRequired('patchSize',@isscalar);
p.addRequired('nRowBipolars',@isscalar);
p.addRequired('nColBipolars',@isscalar);
p.addRequired('rfDiameterMicrons',@isscalar);

p.addParameter('centerNoise',.15,@isscalar); % in units of nBipolars
p.addParameter('baseLineFiringRate',2.2702,@isscalar); % JRG pulled from ON Parasol 2013_08_19_6

vFunc = @(x)(ismatrix(x) || isempty(x));
p.addParameter('ellipseParams',[],vFunc);  % A,B,rho

p.parse(patchSizeMeters,nRowBipolars,nColBipolars,rfDiameterMicrons,varargin{:});

centerNoiseBipolars = p.Results.centerNoise;  
ellipseParams       = p.Results.ellipseParams;  % (Major, Minor, Orientation)
baseLineFiringRate  = p.Results.baseLineFiringRate;

%% sRF Output Details
% After Chichilnisky & Kalmar, 2002
%
% The ellipse intensity maps are defined by this parameterization
%
%   s_center   =   exp(-0.5 * (x-c)*Q * (x-c)') 
%   s_surround = k*exp(-0.5*r*(x-c)*Q*r*(x-c)') 
%
% Here x - c  is a two-dimensional vector that specifies a spatial location.
% Evaluating the exponential, s_center(x-c) indicates the sensitivity at that
% spatial location.  
% x ranges over all positions
% c is a two-dimensional vector that specifies the center of the RF
% Q is a 2 x 2 symmetric positive semi-definite matrix that specifies the
% elliptical Gaussian shape of the RF center, 
% k is a scalar that specifies the relative strength of the surround
% 1/r is a scalar that specifies the relative size of the surround.

% Hard coded for now.  To eliminate.
extent = 10;      % ratio between sampling size and spatial RF standard dev 
r = 0.75;        % radius ratio between center and surround for DoG
k = 1.032 * r;   % scaling of magnitude of surround
  
%% Converting the spatial units from microns to bipolar samples

% Width of the entire patch of cones (which is equal to the patch of
% bipolars) arrives in meters. We convert to microns and build up a
% height/width representation (row,col)
patchSizeMicrons = patchSizeMeters * 1e6;    % Column (width) in um
patchSizeMicronsXY = [patchSizeMicrons*(nRowBipolars/nColBipolars), patchSizeMicrons]; % Row/Col in um

% Determine the number of RGC samples in the hexagonal mosaic
% nRGC: RGC cells, row col. 
nRGC    = floor(patchSizeMicronsXY ./ rfDiameterMicrons); % number of RGCs in h, v direction
nRGC(2) = floor((2/sqrt(3))*nRGC(2));   % Hex packing related

% We calculate the diameter of each bipolar in units of microns.
bipolarDiameterMicrons = (patchSizeMicronsXY(2) / nColBipolars);

% The RGC rfDiameter in terms of the bipolar array
% JRG to explain why there is the 0.5.
rfDiameterBipolars = 0.5 * rfDiameterMicrons / bipolarDiameterMicrons;

% From here the spatial mosaic dimensions are mainly in units of bipolar spacing

% Centers of hexagonally packed receptive fields in bipolar space
rowCenter = (0:2:nRGC(1)-1)*rfDiameterBipolars; % RGC center row coords in nBipolars
colCenter = (sqrt(3)/2 ) *(0:2:nRGC(2)-1)*rfDiameterBipolars; % RGC center col coords in nBipolars
rows = length(rowCenter);
cols = length(colCenter);     % number of RGCs

% This is the sampling range that we use to specify the spatial extent of
% the bipolar cells feeding into one RGC.
pts = -extent*rfDiameterBipolars+1 : extent*rfDiameterBipolars;

%% Create spatial RGC RFs

% Make sure the centers are symmetric around zero.
centerCorrectY = (colCenter(end) - (colCenter(1)))/2; % nBipolars
centerCorrectX = (rowCenter(end) - (rowCenter(1)))/2; % nBipolars

% pre-allocate memory
cellCenterLocations = cell(rows, cols);
spatialRFArray      = cell(rows, cols);
sRFcenter           = cell(rows, cols);
sRFsurround         = cell(rows, cols);
tonicDrive          = cell(rows, cols);

%% Set the tonic drive

% Tonic drive is the bias or DC term for the linear output of the
% GLM. If the tonic drive term is greater than 0, then there is a
% baseline firing rate even when the stimulus input is zero.
% Units of conditional intensity
for ii = 1 : rows
    for jj = 1 : cols
        tonicDrive{ii,jj} = baseLineFiringRate; % from ON Parasol 2013_08_19_6
    end
end

%% Create spatial RFs for each cell

%  The shape of the RF is modeled as an ellipse. 

% Specify the ellipse parameters for each cell

% Jitter the center positions of each cell. 
%
% N.B. Sometimes this introduces a little flip in position.  We could
% eliminate that by using rand() instead of randn().  But we are
% interested.
centerNoiseBipolarsRow = (centerNoiseBipolars*randn(rows,cols))*rfDiameterBipolars;
centerNoiseBipolarsCol = (centerNoiseBipolars*randn(rows,cols))*rfDiameterBipolars;

% These are the ellipse shape parameters (not centered)
ellipseParams = ellipseGen(rows,cols,p.Unmatched,'ellipseParams',ellipseParams);

% The retured ellipse parameters
Qout = cell(rows,cols); 

hexOffset = 0.5 * rfDiameterBipolars;

for ii = 1 : rows
    for jj = 1 : cols
        
        % Compute 2D spatial RF
        
        % Specify RGC centers in bipolar samples.
        % Offset even columns to set the hexagonal packing.
        % Add some jitter to the center positions.  
        % columnOffset = (mod(jj, 2) - 0.5) * rfDiameterBipolars;
        % ic = centerX(ii) + centerNoiseBipolarsRow(ii,jj) - columnOffset;
        % jc = centerY(jj) + centerNoiseBipolarsCol(ii,jj); 
        
        ic = rowCenter(ii) + centerNoiseBipolarsRow(ii,jj);
        jc = colCenter(jj) + centerNoiseBipolarsCol(ii,jj); 
        if mod(jj,2), ic = ic - hexOffset;   % Odd
        else,         ic = ic + hexOffset;   % Even
        end
        
        % JRG - Should we keep the sizes all about the same, or do we allow
        % a lot of variation in the rf sizes.  The norm() operation here
        % forces the overall size of the ellipses to be about 1 bipolar
        % size because the first two dimensions are the major and minor
        % axes.
        if size(ellipseParams,1) == 1
            % User passed in a vector
            ellipseParameters = [ellipseParams(1:2)./norm(ellipseParams(1:2)), ellipseParams(3)];
        else
            % Produced by ellipseGen
            ellipseParameters = [ellipseParams{ii,jj}(1:2)./norm(ellipseParams{ii,jj}(1:2)), ellipseParams{ii,jj}(3)];
        end
        
        % Makes the 2x2 positive definite quadratic form (matrix)
        Qe = ellipseQuadratic(ellipseParameters);
        
        Q = (.125/rfDiameterBipolars^2)*Qe./norm(Qe(:));
                
        % Calculate (x,y) values for input to DoG function in an efficient way
        [i2, j2] = meshgrid(pts, pts); % nBipolars
        i = i2(:); j = j2(:);          % nBipolars
        IJ = [i j];
        
        % Scale by the r and Q
        QIJ  = Q*IJ'; 
        rQIJ = r*Q*IJ';       % unitless
        %  icrm = repmat([ic jc],length(i),1);
        
        % (-0.5*(x-c)*Q*(x-c)'): unitless
        p1 = prod([IJ(:,1) QIJ(1,:)'],2) + prod([IJ(:,2) QIJ(2,:)'],2);
        p2 = prod([IJ(:,1) rQIJ(1,:)'],2)+ prod([IJ(:,2) rQIJ(2,:)'],2);
        
        % DoG calculation
        % conditional intensity, related by Poisson firing to spikes/sec
        so_center = reshape(exp(-0.5*p1), size(i2));
        so_surround = reshape(k*exp(-0.5*p2), size(i2));
        so = so_center - so_surround;
        
        % Store calculated parameters, units of conditional intensity
        cellCenterLocations{ii,jj} = (patchSizeMicronsXY(2) / nColBipolars)*([ic jc] - [centerCorrectX centerCorrectY]); % nBipolars
        spatialRFArray{ii,jj} = so;
        sRFcenter{ii,jj} = so_center;
        sRFsurround{ii,jj} = so_surround;
        
        % Do some calculations to make plots where RFs are filled in
        % Measure magnitude at 1 SD from center
        if ii == 1 && jj == 1
            xv = [1 0];   % rand(1,2);
            xvn = rfDiameterBipolars * xv./norm(xv);
            x1 = xvn(1); y1 = xvn(2);
            magnitude1STD = exp(-0.5*[x1 y1]*Q*[x1; y1])- k*exp(-0.5*[x1 y1]*r*Q*[x1; y1]);
            [maxv,maxr] = max(so_center(:)-so_surround(:)); [mr,mc] = ind2sub(size(so_center),maxr);
            rii = mr; cii = mc; im = 1;
            while (so_center(mr,cii)-so_surround(mr,cii)) > magnitude1STD; im = im+1; cii = mc-1+im; end; [rfDiameterBipolars (cii-mc-1)]
        end
        Qout{ii,jj} = ellipseParameters;
    end
end

end