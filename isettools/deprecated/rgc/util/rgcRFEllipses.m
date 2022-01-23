function [sRFcenter, sRFsurround, cellLocation, ellipseMatrix] = ...
    rgcRFEllipses(cellLocation, spread, varargin)
% build array of elliptical receptive fields
%
% Syntax:
%   [sRFcenter, sRFsurround, cellLocation, ellipseMatrix] = ...
%       rgcRFEllipses(cellLocation, rfDiameter, varargin)
%
% Description:
%    The ellipse intensity maps are defined by this parameterization
%       s_center = exp(-0.5 * (x - c) * Q * (x - c)')
%       s_surround = k * exp(-0.5 * r * (x - c) * Q * r * (x - c)')
%
%    Here x - c  is a 2D vector that specifies a spatial location.
%    Evaluating the exponential, s_center(x - c) indicates the sensitivity
%    at that spatial location.
%    x ranges over all positions
%    c is a two-dimensional vector that specifies the center of the RF
%    Q is a 2 x 2 symmetric positive semi-definite matrix that specifies
%    the elliptical Gaussian shape of the RF center, 
%    k is a scalar that specifies the relative strength of the surround
%    1 / r is a scalar that specifies the relative size of the surround.
%
% Inputs:
%    cellLocation            - Matrix. A matrix denoting cell locations, in
%                              units of bipolar samples.
%    spread                  - Numeric. A single standard deviation in
%                              units of bipolar samples.
%
% Outputs:
%    sRFcenter               - Cell. A nRows by nCols cell array containing
%                              the sRFcenter. 
%    sRFsurround             - Cell. A nRows by nCols cell array containing
%                              the sRFsurround.
%    cellLocation            - Matrix. The modified cellLocation matrix.
%    ellipseMatrix           - Matrix. an ellipse generated using the
%                              available ellipse parameters.
%
% Optional key/value pairs:
%    centerSurroundSizeRatio - Numeric. The center surround size ratio.
%                              Default is sqrt(0.75)
%    centerSurroundAmpRatio  - Numeric. The center surround amplitude
%                              ratio. Default is 0.6 * 0.774.
%    centerNoise             - Numeric. The center noise. Default is 0.15.
%    ellipseParams           - Matrix. A matrix of [A, B, rho]. Default is
%                              an empty matrix, [].
%
% References
%    * spatial RF  Details - After Chichilnisky & Kalmar, 2002
%

% History:
%    XX/XX/17  JRG/BW  ISETBIO Team, 2017
%    05/31/19  JNM     Documentation pass

%%
p = inputParser;

% cellLocation - in units of bipolar samples
p.addRequired('cellLocation', @isnumeric);
p.addRequired('spread', @isscalar);  % 1 std. in units of bipolar samples

% JRG has a theoretical argument about this and we will fill it in. But for
% now, these are hard coded numerical defaults that satisfy the
% descriptions in the paper by Chichilnisky and Kalmar, 2002 (around page
% 2738).
%
% r = sqrt(0.75);     % radius ratio between center and surround for DoG
% k = 1.032 * r ^ 2;  % scaling of magnitude of surround
%
p.addParameter('centerSurroundSizeRatio', sqrt(0.75), @isscalar);
p.addParameter('centerSurroundAmpRatio', .6 * 0.774, @isscalar);
p.addParameter('centerNoise', .15, @isscalar);  % bipolar sample units
% ellipseParams - A, B, rho
p.addParameter('ellipseParams', [], @(x)(ismatrix(x) || isempty(x)));
p.parse(cellLocation, spread, varargin{:});

centerNoiseBipolars = p.Results.centerNoise;  % In bipolar sample units
ellipseParams = p.Results.ellipseParams;      % (Major, Minor, Orientation)
% r - The radius ratio between center and surround for DoG
r = p.Results.centerSurroundSizeRatio;
k = p.Results.centerSurroundAmpRatio;         % surround magnitude scaling

%% Allocate space
[nRows, nCols, ~] = size(cellLocation);

% Each spatial RF can differ a bit.  So, we make them a cell array
sRFcenter = cell(nRows, nCols);
sRFsurround = cell(nRows, nCols);

%% Create the centers of the RGC sample positions

% Specify the spatial extent of the bipolar samples feeding into one RGC.
% We make them extend +/- three standard deviations (spread is 1 std dev).
pts = (-2*spread):(2*spread);

%% Jitter the center positions of each cell.
% N.B. Sometimes this introduces a little flip in position. We could
% eliminate that by using rand() instead of randn().
cellLocation = cellLocation + ...
    centerNoiseBipolars * randn(nRows, nCols, 2) * spread;

% These are the ellipse shape parameters (not centered)
ellipseMatrix = ellipseGen(nRows, nCols, p.Unmatched, ...
    'ellipseParams', ellipseParams);

% To build the hex mosaic, we need this value
hexOffset = 0.5 * spread;

% Compute 2D spatial RF (ellipses)
for rr = 1:nRows      % Row index
    for cc = 1:nCols  % Col index
        % Specify RGC centers in bipolar sample space.
        % First, add some jitter to the center positions.
        thisRowCenter = cellLocation(rr, cc, 1);
        thisColCenter = cellLocation(rr, cc, 2); 

        % Offset every other column to create a hexagonal packing.
        if mod(cc, 2)  % Odd case
            thisRowCenter = thisRowCenter - 1 * hexOffset / 1;
        else           % Even case
            thisRowCenter = thisRowCenter;  % + hexOffset / 2;
        end
        % Save the new cell center location in bipolar samples
        cellLocation(rr, cc, :) = [thisRowCenter, thisColCenter];

        % Makes the 2x2 positive definite quadratic form (matrix)
        % In order to keep the same area under the DoG surface, need to
        % normalize the diagonal.
        % Q = (1 / (.5 * spread) ^ 2) * diag(ellipseMatrix{rr, cc}(1:2));
        Q = (1 / ((.25 * spread) ^ 2)) * diag(ellipseMatrix{rr, cc}(1:2));

        % For the DoG, we need to do the rotation matrix separately from Q,
        % otherwise the DoG height and width change for the same magnitude
        % params with a different angle param.
        R = [cosd(ellipseMatrix{rr, cc}(3)), ...
            -sind(ellipseMatrix{rr, cc}(3));...
            sind(ellipseMatrix{rr, cc}(3)), ...
            cosd(ellipseMatrix{rr, cc}(3))];

        % Calculate the (x, y) values for input to DoG function in an
        % efficient way.
        [X, Y] = meshgrid(pts, pts); % units of bipolar samples

        % Apply rotation matrix to the receptive field
        XY = (R * [X(:) Y(:)]')';

        % The calculation we are working up to is this:
        %              (-0.5 * x * Q * x')
        % Which is used for the Difference of Gaussian RF shape
        Qmatr = Q * XY';
        rQmatr = r ^ 2 * Q * XY';  % unitless

        % The prod is faster than other ways
        QXY = prod([XY(:, 1)  Qmatr(1, :)'], 2) + ...
            prod([XY(:, 2)  Qmatr(2, :)'], 2);
        RQXY = prod([XY(:, 1) rQmatr(1, :)'], 2) + ...
            prod([XY(:, 2) rQmatr(2, :)'], 2);

        % DoG calculation
        % conditional intensity, related by Poisson firing to spikes/sec
        so_center = reshape(exp(-0.5 * QXY), size(X));
        so_surround = reshape(k * exp(-0.5 * RQXY), size(X));

        so_surf = so_center-so_surround;
        normFactor = 2 / sum(so_surf(:));
        % spatialRFArray{ii, jj} = so;
        % Store calculated parameters, units of conditional intensity
        sRFcenter{rr, cc} = normFactor * so_center;
        sRFsurround{rr, cc} = normFactor * so_surround;
    end
end

end
