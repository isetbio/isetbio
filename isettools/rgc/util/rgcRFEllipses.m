function [sRFcenter, sRFsurround, cellLocation, ellipseMatrix] = rgcRFEllipses(nRGC,rfDiameter,varargin)
% RGCRFELLIPSES - build array of elliptical receptive fields
%
%   [sRFcenter, sRFsurround, cellLocation, ellipseMatrix] = ...
%       rgcRFEllipses(nRGC,rfDiameter,varargin)
%
% Inputs
%
% Optional inputs
%
%
% spatial RF  Details
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
%
% JRG/BW ISETBIO Team, 2017

%% Programming:  See end for a test

%%
p = inputParser;

p.addRequired('nRGC',@isvector); % in units of nBipolars
p.addRequired('rfDiameter',@isscalar); % in units of nBipolars

p.addParameter('centerNoise',.15,@isscalar); % in units of nBipolars
p.addParameter('baseLineFiringRate',2.2702,@isscalar); % JRG pulled from ON Parasol 2013_08_19_6

vFunc = @(x)(ismatrix(x) || isempty(x));
p.addParameter('ellipseParams',[],vFunc);  % A,B,rho
p.parse(nRGC,rfDiameter,varargin{:});

centerNoiseBipolars = p.Results.centerNoise;    % In bipolar sample units
ellipseParams       = p.Results.ellipseParams;  % (Major, Minor, Orientation)

%%
% Hard coded for now.  To eliminate.
extent = 2.5;        % ratio between sampling size and spatial RF standard dev 
r = sqrt(0.75);      % radius ratio between center and surround for DoG
k = 1.032 * r^2;     % scaling of magnitude of surround

%%
% The potentially jittered cell centers 
cellLocation = zeros(nRGC(1), nRGC(2),2);

% Each spatial RF can differ a bit.  So, we make them a cell array
sRFcenter           = cell(nRGC(1),nRGC(2));
sRFsurround         = cell(nRGC(1),nRGC(2));

nRows = nRGC(1); nCols = nRGC(2);

%%
% Initial, noise free sample positions of RGC receptive field centers in
% bipolar space.  The centers are spaced by the one rfDiameter and the
% coordinate frame has a (0,0) in the middle of the bipolar plane.  No
% noise at the moment.  Added later.
rowCenters = (0:nRows-1)*rfDiameter;               
colCenters = (sqrt(3)/2 ) *(0:nCols-1)*rfDiameter;
rowCenters = rowCenters - mean(rowCenters);
colCenters = colCenters - mean(colCenters);

% This is the sampling range that we use to specify the spatial extent of
% the bipolar cells feeding into one RGC.  This should be bigger than the
% largest rgc RF.  So really, extent should be chosen based on the sizer of
% the RGC RFs, not fixed as it is here (BW).
pts = -extent*rfDiameter:extent*rfDiameter;

%%
% Jitter the center positions of each cell.
% N.B. Sometimes this introduces a little flip in position.  We could
% eliminate that by using rand() instead of randn().
centerNoiseRows = (centerNoiseBipolars*randn(nRows,nCols))*rfDiameter;
centerNoiseCols = (centerNoiseBipolars*randn(nRows,nCols))*rfDiameter;
% vcNewGraphWin; plot(centerNoiseBipolarsCol(:),centerNoiseBipolarsRow(:),'.')

% These are the ellipse shape parameters (not centered)
ellipseMatrix = ellipseGen(nRows,nCols,p.Unmatched,'ellipseParams',ellipseParams);

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
        % ellipseP = [ellipseParams{rr,cc}(1:2), ellipseParams{rr,cc}(3)];
        
        % Makes the 2x2 positive definite quadratic form (matrix)
        % In order to keep the same area under the DoG surface, need to
        % normalize the diagonal.
        Q =  (1/(.5*rfDiameter)^2)*diag(ellipseMatrix{rr,cc}(1:2));
        
        % For the DoG, we need to do the rotation matrix separately from Q,
        % otherwise the DoG height and width change for the same magnitude
        % params with a different angle param.
        R = [cosd(ellipseMatrix{rr,cc}(3)) -sind(ellipseMatrix{rr,cc}(3));...
            sind(ellipseMatrix{rr,cc}(3))   cosd(ellipseMatrix{rr,cc}(3))];
        
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
        
        % Save the cell center location in bipolar samples
        cellLocation(rr,cc,:) = [thisRowCenter, thisColCenter];
        
        % spatialRFArray{ii,jj} = so;
        % Store calculated parameters, units of conditional intensity
        sRFcenter{rr,cc}      = so_center;
        sRFsurround{rr,cc}    = so_surround;
        
    end
end

end


% Since we create the plot of the RF mosaic as an ellipse, but the
% actual RFs as DoGs, we need to check that the DoG magnitude at
% rfDiameter is actually 1 std. We do that here by first
% calculating what the 1 std magnitude of the DoG should be, and
% then comparing it to what our actual RF has. They should match
% within 1 (units of bipolar samples). If there is high variance
% in the shapes of RFs, then individual RFs might not match, but
% they should on average.
% To find the radius at which our 1 std magnitude occurs in our DoG
% RF, we find the coordinates of the max values of the RF. Then we
% move one bipolar sample away from this and check if the value has
% decreased below the specified 1 std magnitude; if not, then move
% another bipolar sample away, and do this again until the weight
% on the bipolar sample is below the 1 std magnitude. Then we
% compare our observed value of radius in bipolar samples to the
% pre-specified radius.

%         if rr == 1 && cc == 1
%             % Find the RF weight at a distance of the diameter
%             xv = [1 0];
%             xvn = rfDiameter * xv./norm(xv);
%             x1 = xvn(1); y1 = xvn(2);
%             % Calculate value at 1 std
%             magnitude1STD = exp(-0.5*[x1 y1]*Q*[x1; y1])- k*exp(-0.5*[x1 y1]*r^2*Q*[x1; y1]);
%             % Find coordinates at 1 std
%             [~,maxr] = max(so_center(:)-so_surround(:)); [mr,mc] = ind2sub(size(so_center),maxr);
%             cii = mc; im = 1;
%             % Move one sample away until we find the 1 std sample
%             while (cii < size(so_center,2)) && ((so_center(mr,cii)-so_surround(mr,cii)) > magnitude1STD)
%                 im = im+1;
%                 cii = mc-1+im;
%             end
%             % [rfDiameter (cii-mc-1)]
%             % if abs(rfDiameter - (cii - mc - 1)) > 1; display('RF mismatch');
%         end
%