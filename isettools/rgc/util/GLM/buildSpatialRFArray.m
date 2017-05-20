function [sRFcenter, sRFsurround, cellCenterLocations, tonicDrive, ellipseParams] = ...
    buildSpatialRFArray(patchSizeMeters, nRowBipolars, nColBipolars, rfDiameterMicrons, varargin)
% BUILDSPATIALARRAY - Create RF center and surround for RGC mosaic
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
extent = 2.5;     % ratio between sampling size and spatial RF standard dev 
r = sqrt(0.75);        % radius ratio between center and surround for DoG
k = 1.032 * r^2;   % scaling of magnitude of surround
  
%% Converting the spatial units from microns to bipolar samples

% p[atchSizeMeters is the width of the patch of cones sampled by the
% bipolars. This arrives in meters, we convert row/col in microns
patchRowColMicrons = [patchSizeMeters*(nRowBipolars/nColBipolars), patchSizeMeters]*1e6; % Row/Col in um

% This is the diameter of the cone mosaic sampled by each bipolar in units
% of microns.
bipolarDiameterMicrons = (patchRowColMicrons(2) / nColBipolars);

%% The spatial mosaic dimensions are in units of bipolar spacing

% The RGC rfDiameter in terms of bipolar array samples.  So, one step is
% one step in the bipolar array. From here on out, the spatial coordinate
% system is with respect to the bipolar array.
rfDiameter = rfDiameterMicrons / bipolarDiameterMicrons;

% Determine the number of RGC samples in the hexagonal mosaic
% nRGC: RGC cells, row col. 
nRGC    = floor(patchRowColMicrons ./ rfDiameterMicrons); % number of RGCs in h, v direction

% Adjust the scale factor to account for the hexagonal packing of the RGC
% mosaic
nRGC(2) = floor((2/sqrt(3))*nRGC(2));

% Initial, noise free sample positions of RGC receptive field centers in
% bipolar space.  The centers are spaced by the one rfDiameter and the
% coordinate frame has a (0,0) in the middle of the bipolar plane.  No
% noise at the moment.  Added later.
rowCenters = (0:nRGC(1)-1)*rfDiameter;               
colCenters = (sqrt(3)/2 ) *(0:nRGC(2)-1)*rfDiameter;
rowCenters = rowCenters - mean(rowCenters);
colCenters = colCenters - mean(colCenters);

% Count the number of rows and columns of the RGCs
nRows = length(rowCenters);
nCols = length(colCenters);

% This is the sampling range that we use to specify the spatial extent of
% the bipolar cells feeding into one RGC.  This should be bigger than the
% largest rgc RF.  So really, extent should be chosen based on the sizer of
% the RGC RFs, not fixed as it is here (BW).
pts = -extent*rfDiameter:extent*rfDiameter;

%% Create spatial RGC RFs

% Make sure the centers are symmetric around zero.
% centerCorrectY = (colCenters(end) - (colCenters(1)))/2; % nBipolars
% centerCorrectX = (rowCenters(end) - (rowCenters(1)))/2; % nBipolars

% pre-allocate memory
cellCenterLocations = cell(nRows, nCols);
sRFcenter           = cell(nRows, nCols);
sRFsurround         = cell(nRows, nCols);
tonicDrive          = cell(nRows, nCols);

%% Set the tonic drive

% Tonic drive is the bias or DC term for the linear output of the
% GLM. If the tonic drive term is greater than 0, then there is a
% baseline firing rate even when the stimulus input is zero.
% Units of conditional intensity
for rr = 1 : nRows
    for cc = 1 : nCols
        tonicDrive{rr,cc} = baseLineFiringRate; % from ON Parasol 2013_08_19_6
    end
end

%% Create spatial RFs for each cell

%  The shape of the RF is modeled as an ellipse. 

% Specify the ellipse parameters for each cell

% Jitter the center positions of each cell. 
% N.B. Sometimes this introduces a little flip in position.  We could
% eliminate that by using rand() instead of randn().
centerNoiseRows = (centerNoiseBipolars*randn(nRows,nCols))*rfDiameter;
centerNoiseCols = (centerNoiseBipolars*randn(nRows,nCols))*rfDiameter;
% vcNewGraphWin; plot(centerNoiseBipolarsCol(:),centerNoiseBipolarsRow(:),'.')

% These are the ellipse shape parameters (not centered)
ellipseParams = ellipseGen(nRows,nCols,p.Unmatched,'ellipseParams',ellipseParams);

% The retured ellipse parameters
% Qout = cell(rows,cols); 

% To build the hex mosaic, we need this value
hexOffset = 0.5 * rfDiameter;

for rr = 1:nRows         % Row index
    for cc = 1:nCols     % Col index
        
        % Compute 2D spatial RF
        
        % Specify RGC centers in bipolar sample space.
        % First, add some jitter to the center positions.  
        thisRowCenter = rowCenters(rr) + centerNoiseRows(rr,cc);
        thisColCenter = colCenters(cc) + centerNoiseCols(rr,cc);
        
        % Then offset columns to set the hexagonal packing.
        if mod(cc,2), thisRowCenter = thisRowCenter - hexOffset;   % Odd
        % else,         thisRowCenter = thisRowCenter + hexOffset;   % Even
        end        
        
        % Without the normalization
        ellipseP = [ellipseParams{rr,cc}(1:2), ellipseParams{rr,cc}(3)];

        % Makes the 2x2 positive definite quadratic form (matrix)        
        % In order to keep the same area under the DoG surface, need to
        % normalize the diagonal.
        Q =  (1/(.5*rfDiameter)^2)*diag(ellipseParams{rr,cc}(1:2));                 
                
        % For the DoG, we need to do the rotation matrix separately from Q,
        % otherwise the DoG height and width change for the same magnitude
        % params with a different angle param.
        R = [cosd(ellipseParams{rr,cc}(3)) -sind(ellipseParams{rr,cc}(3));...
            sind(ellipseParams{rr,cc}(3))   cosd(ellipseParams{rr,cc}(3))];
        
        % Calculate (x,y) values for input to DoG function in an efficient way
        [X, Y] = meshgrid(pts, pts); % nBipolars
        % XY = [X(:) Y(:)];       
        % Apply rotation matrix; need to do it here to coordinates so that
        % the magnitude of the DoG is not incorrectly scaled.
        XY = (R*[X(:) Y(:)]')';        
        
        % % % This is very slow
        % Scale by the r and Q
        % QXY  = diag(XY * Q * XY'); %         
        % % Surround
        % RQXY = r^2*QXY;       % unitless        
        
        % % % % This does the same thing but much faster for big RFs
        Qmatr  = Q*XY'; 
        rQmatr = r^2*Q*XY';       % unitless        
        % (-0.5*(x-c)*Q*(x-c)'): unitless
        QXY =  prod([XY(:,1)  Qmatr(1,:)'],2) + prod([XY(:,2)  Qmatr(2,:)'],2);
        RQXY = prod([XY(:,1) rQmatr(1,:)'],2) + prod([XY(:,2) rQmatr(2,:)'],2);
        
        % DoG calculation
        % conditional intensity, related by Poisson firing to spikes/sec
        so_center   = reshape(exp(-0.5*QXY),    size(X));
        so_surround = reshape(k*exp(-0.5*RQXY), size(X));        
        
        % Save the cell center location
        cellCenterLocations{rr,cc} = [thisRowCenter, thisColCenter];
        
        % spatialRFArray{ii,jj} = so;
        % Store calculated parameters, units of conditional intensity
        sRFcenter{rr,cc}      = so_center;
        sRFsurround{rr,cc}    = so_surround;
        
        % Since we create the plot of the RF mosaic as an ellipse, but the
        % actual RFs as DoGs, we need to check that the DoG magntiude at
        % rfDiameter is actually 1 std. We do that here by first
        % calculating what the 1 std magnitude of the DoG should be, and
        % then comparing it to what our actual RF has. They should match
        % within 1 (units of bipolar samples). If there is high variance
        % in the shapes of RFs, then individual RFs might not match, but
        % they should on average.
        if rr == 1 && cc == 1
            xv = [1 0];   % rand(1,2);
            xvn = rfDiameter * xv./norm(xv);
            x1 = xvn(1); y1 = xvn(2);
            magnitude1STD = exp(-0.5*[x1 y1]*Q*[x1; y1])- k*exp(-0.5*[x1 y1]*r^2*Q*[x1; y1]);
            [maxv,maxr] = max(so_center(:)-so_surround(:)); [mr,mc] = ind2sub(size(so_center),maxr);
            rii = mr; cii = mc; im = 1;
            while (so_center(mr,cii)-so_surround(mr,cii)) > magnitude1STD; im = im+1; cii = mc-1+im; end; [rfDiameter (cii-mc-1)]
        end
        
    end
end

end