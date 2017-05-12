function [sRFcenter, sRFsurround, cellCenterLocations, tonicDrive, Qout] = ...
    buildSpatialRFArray(patchSizeMeters, nRowBipolars, nColBipolars, rfDiameterMicrons, varargin)
% Builds the spatial RF center and surround for the cells in a mosaic
% 
% The spatial RFs are generated according to the size of the pixel, cone or
% bipolar mosaic, their spacing (in microns) and the diameter of the RGC RF
% as determined by the TEE of the retinal patch.
% 
%  [sRFcenter, sRFsurround, rfDiaMagnitude, cellCenterLocations] = 
%      buildSpatialRFArray(spacing, row, col, rfDiameter)
% 
% Inputs: 
%       spacing - Center to center of the RF in microns
%       row     - Number of input samples
%       col     - Number of input samples
%       rfDiameter - receptive field of 1 std in microns
%   
% Outputs: 
%       spatialRFcenter cell array, 
%       spatialRFsurround cell array,
%       cellCenterLocations cell array.
%       tonicDrive cell array
%           Tonic drive is the bias or DC term for the linear output of the GLM. 
% 
% Example:
%   Build spatial RFs of the RGCs in this mosaic
%
% [obj.sRFcenter, obj.sRFsurround, obj.rfDiaMagnitude, obj.cellLocation] = ...
%     buildSpatialRFArray(innerRetina.spacing, innerRetina.row, innerRetina.col, obj.rfDiameter);
% 
% 9/2015 JRG (c) isetbio
% 7/2016 JRG updated

%% Manage parameters
p = inputParser;
p.addRequired('patchSize',@isscalar);
p.addRequired('nRowBipolars',@isscalar);
p.addRequired('nColBipolars',@isscalar);
p.addRequired('rfDiameterMicrons',@isscalar);
p.addParameter('centerNoise',.15,@isscalar); % in units of nBipolars
p.addParameter('baseLineFiringRate',2.2702,@isscalar); % JRG pulled from ON Parasol 2013_08_19_6
p.KeepUnmatched = true;
vFunc = @(x)(ismatrix(x) || isempty(x));
p.addParameter('ellipseParams',[],vFunc);  % A,B,rho

p.parse(patchSizeMeters,nRowBipolars,nColBipolars,rfDiameterMicrons,varargin{:});

centerNoiseBipolars = p.Results.centerNoise;  
ellipseParams       = p.Results.ellipseParams;  % (Major, Minor, Orientation)
baseLineFiringRate  = p.Results.baseLineFiringRate;

%% sRF Output Details
% s_center   =   exp(-0.5 * (x-c)*Q * (x-c)') 
% s_surround = k*exp(-0.5*r*(x-c)*Q*r*(x-c)') 

% % Chichilnisky & Kalmar, 2002
% Here x is a two-dimensional vector that specifies a spatial location, s(x)
% indicates the sensitivity at that spatial location, c is a
%  two-dimensional vector that specifies the midpoint of the RF, Q is a 2 x
% 2 symmetric positive semi-definite matrix that specifies the elliptical
% Gaussian shape of the RF center, k is a scalar that specifies the
% relative strength of the surround, and 1/r is a scalar that specifies the
% relative size of the surround.

% Hard coded for now.  To eliminate.
extent = 5;      % ratio between sampling size and spatial RF standard dev 
r = 0.75;        % radius ratio between center and surround for DoG
k = 1.032 * r;   % scaling of magnitude of surround
  
%%

% I think JRG built this assuming the patch size referred to a region of
% the cone mosaic.  
%
% Size of bipolar patch arrives in meters, which apparently is the width of
% the patch.  We convert to microns and build up a height/width
% representation (row,col)
patchSizeMicrons = patchSizeMeters * 1e6;                % Column (width) in um

% Adjust aspect ratio of RGC array if bipolar mosaic is not square
patchSizeMicronsXY = [patchSizeMicrons*(nRowBipolars/nColBipolars), patchSizeMicrons]; % Row/Col in um

% Determine the number of RGC samples in the mosaic
% patchSizeMicrons: um; rfDiameterMicrons: um/RGC cell; nRGC: RGC cells 
nRGC = floor(patchSizeMicronsXY ./ rfDiameterMicrons); % number of RGCs in h, v direction
nRGC(2) = floor((2/sqrt(3))*nRGC(2));

% The rfDiameter comes here in units of um, and this 
% converts the rf diameter to units of number of bipolars
% rfDiameter in: um/RGC, 
% patchSize: um, inCol: number bipolar cells
% (patchSize(2) / inCol) : um/bipolar cell
% rfDiameter out: number bipolar cells per RGC
rfDiameterBipolars = 1*rfDiameterMicrons / (patchSizeMicronsXY(2) / nColBipolars);

% centers of receptive fields
centerX = (0:2:nRGC(1)-1)*rfDiameterBipolars; % RGC center row coords in nBipolars
centerY = (sqrt(3)/2 ) *(0:2:nRGC(2)-1)*rfDiameterBipolars; % RGC center col coords in nBipolars
rows = length(centerX); cols = length(centerY);     % number of RGCs

% number bipolar cells out to the extent of the spatial RF
pts = -extent*rfDiameterBipolars+1 : extent*rfDiameterBipolars;

%% Create spatial RFs
% offset for zero centering in number of bipolars
centerCorrectY = (centerY(end) - (centerY(1)))/2; % nBipolars
centerCorrectX = (centerX(end) - (centerX(1)))/2; % nBipolars

% pre-allocate space
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
%  Effectively creates ellipse quadratic, 
%  Center positions and general sRF() receptive field weightings functions

% Compute the variation in the center position for every cell
% This could flip the position occasionally.  We were using
% rand(row,col)*2 - 1 to guarantee no flips.
centerNoiseBipolarsRow = centerNoiseBipolars*rfDiameterBipolars*(randn(rows,cols));
centerNoiseBipolarsCol = centerNoiseBipolars*rfDiameterBipolars*(randn(rows,cols));

Qout = cell(rows,cols); 

ellipseParams = ellipseGen(rows,cols,p.Unmatched,'ellipseParams',ellipseParams);
 
for ii = 1 : rows
    for jj = 1 : cols
        
        % Compute 2D spatial RF
        
        % Specify centers in um, offset even rows for hexagonal packing
        % Add some jitter to the center positions.  The size of the jitter
        % is the parameter centerNoiseBipolars
        ic = centerX(ii) + centerNoiseBipolarsRow(ii,jj) - (mod(jj, 2) - 0.5) * rfDiameterBipolars;
        jc = centerY(jj) + centerNoiseBipolarsCol(ii,jj); 
   
        if size(ellipseParams,1)==1
            ellipseParameters = [ellipseParams(1:2)./norm(ellipseParams(1:2)), ellipseParams(3)];
        else
            ellipseParameters = [ellipseParams{ii,jj}(1:2)./norm(ellipseParams{ii,jj}(1:2)), ellipseParams{ii,jj}(3)];
        end
        ellipseParameters
        Qe = ellipseQuadratic(ellipseParameters); 
        Q = (.125/rfDiameterBipolars^2)*Qe./norm(Qe(:));
                
        % Calculate (x,y) values for input to DoG function in an efficient way
        [i2, j2] = meshgrid(pts, pts); % nBipolars
        i = i2(:); j = j2(:);          % nBipolars
        IJ = [i j];
        
        % Scale by the r and Q
        QIJ = Q*IJ'; 
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
            xv = [1 0];%rand(1,2);
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