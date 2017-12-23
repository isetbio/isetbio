function [zzgrid, xvec] = ffndgrid(x, f, delta, limits, aver)
% Fast 'n' Furious N-D data gridding.
%
% Syntax:
%   [fgrid, xvec] = ffndgrid(x, f, delta, limits, aver);
%
% Description:
%    Grids unevenly spaced data in vector f into a matrix grid.
%
% Inputs:
%    x      - [x1 x2 ... xD] coordinate matrix for unevenly spaced data, f.
%             Size N x D. 
%    f      - f(x), vector of function values length N.
%    delta  - [dx1 + i * pad, dx2, ..., dxD] or [-Nx1 + i * pad, -Nx2, ...
%             NxD], where dx1 to dxD and Nx1 to NxD defines the stepsize of
%             grids and number of bins, respectively, in the D dimesional
%             space. Empty gridpoints are padded with IMAG(delta(1)) = pad.
%             (default -75)
%   limits  - [xl1 xu1 xl2 ...xuN fl fu n0], contain the limits in the
%             D-dimensional x1-x2...xN-f-plane of where to grid. The last
%             parameter, n0, removes outliers from the data set by ignoring
%             grid points with n0 or less observations. When n0 is negative
%             it is treated as the percentage of the total number of data
%             points. (default [min(x1), max(x1), min(x2), .... max(xN),
%             min(f), max(f), 0])
%    aver   - (Optional) Determine how to calculate the values of fgrid.
%             Possible values of 0 and 1. Default is 1. Each value of fgrid
%             will be calculated by performing the specified operation on
%             all of the points within the cell. The values correlate to
%             the following operations:
%                 0 - Sum
%                 1 - Average (Default)
%
% Outputs:
%    zzgrid - Matrix of gridded data.
%    xvec   - cellarray of gridvectors. This will be xl1:dx1:xu1 or
%             linspace(xl1, xu1, Nx1) depending on delta.
%
% Notes:
%    * [Note:  BW External function.  We aren't worrying about it.  Probably it
%    should be in external.]
%    * [Note: XXX - The vector limits can be padded with NaNs if only
%      certain limits are desired, e g if xl1 and fl are wanted:
%            ffndgrid(x, f, [], [.5 nan nan nan 45])
%      If no output arguments are given, FFGRID will plot the gridded
%      function with the prescribed axes using pColor.]
%    * [Note: XXX - Tested on: MatLab 4.2, 5.0, 5.1, 5.2 and 5.3.]
%    * [Note: JNM - ~Line 214 - there is an 'if 1' statement with an else
%      case, why is this unreachable code included, or should there have
%      been additional calculation done in the if statement rather than a
%      universal true?]
%
% See Also:
%    griddata
%

% History:
%    07/28/97  OB   Oyvind.Breivik@gfi.uib.no.
%                   Oyvind Breivik, Department of Geophysics,
%                   University of Bergen, Norway
%    06/12/98  PAB  ???
%    10/05/98  PAB  secret option: aver. Optionally do not take average of
%                   values for each point
%    02/08/01  PAB  Made it general for D dimensions + changed name from
%                   ffgrid to ffndgrid.
%                   Added nargchk + examples.
%                   Updated help header to wafo-style.
%                   Moved dx and dy into delta = [dx, dy].
%                   Removed call to bin
%    12/01/17  JNM  Formatting

% Examples:
%{
    N = 500;
    D = 2;
    sz = [N, D];
    x = randn(sz);
    z = ones(sz(1), 1); 
    [nc, xv] = ffndgrid(x, z, -15, [], 0);  % Histogram
    pcolor(xv{:}, nc)     %
    [XV, YV] = meshgrid(xv{:});
    text(XV(:), YV(:), int2str(nc(:)))
    dx = [diff(xv{1}(1:2)) diff(xv{2}(1:2))];
    contourf(xv{:}, nc / (N * prod(dx))) % 2-D probability density plot.
    colorbar
    colormap jet
%}

narginchk(2, 5)		

[r, ~] = size(x);
% Make sure x is a column vector.
if r == 1
    x = x(:);
end

[N, D] = size(x);
f = f(:);
if length(f) == 1
    f = f(ones(N, 1), :);
elseif length(f) ~= N
    error('The length of f must equal size(x, 1)!')
end

% Default values
dx = repmat(-75, 1, D);
pad = 0; 
xyz = [zeros(1, 2 * D) min(f), max(f), 0];
xyz(1:2:2 * D) = min(x, [], 1);
xyz(2:2:2 * D) = max(x, [], 1);


% take average of values for each point (default)
if (nargin < 5) || isempty(aver), aver = 1; end
if (nargin >= 4) && ~isempty(limits)
    nlim = length(limits);
    ind = find(~isnan(limits(1:min(7, nlim))));
    if any(ind), xyz(ind) = limits(ind); end
end
if nargin >= 3 && ~isempty(delta)
    Nd = length(delta);
    delta = delta(1:min(Nd, D));
    if Nd == 1, delta = repmat(delta, 1, D); end
    ind = find(~(isnan(delta) | delta == 0));
    if any(ind)
        dx(ind) = real(delta(ind)); 
        pad = imag(delta(1));
    end
end

xL = xyz(1:2:2 * D);
xU = xyz(2:2:2 * D);

fL = xyz(end - 2);
fU = xyz(end - 1);
n0 = xyz(end);

ind = find(dx < 0);
if any(ind)
    if any(dx(ind) ~= round(dx(ind)))
        error('Some of Nx1, ...NxD in delta are not an integer!'), 
    end
    dx(ind) = (xU(ind) - xL(ind)) ./ (abs(dx(ind)) - 1);
end

% bin data in D-dimensional-space
binx = round((x - xL(ones(N, 1), :)) ./ dx(ones(N, 1), :)) + 1;

fgsiz = ones(1, max(D, 2));
xvec = cell(1, D);
for iy=1:D
    xvec{iy} = xL(iy):dx(iy):xU(iy);
    fgsiz(iy) = length(xvec{iy});
end
if D > 1
    in = all((binx >= 1) & (binx <= fgsiz(ones(N, 1), :)), 2) ...
        & (fL <= f) & (f <= fU);
else
    in = (binx >= 1) & (binx <= fgsiz(1)) & (fL <= f) & (f <= fU);
end
binx = binx(in, :); 
f = f(in);
N = length(binx); % how many datapoints are left now?
Nf = prod(fgsiz);

if D>1
    fact = [1 cumprod(fgsiz(1:D - 1))];
    % linear index to fgrid
    binx = sum((binx - 1) .* fact(ones(N, 1), :), 2) + 1;
end

% Fast gridding
fgrid = sparse(binx, 1, f, Nf, 1);% z-coordinate

if n0 ~= 0 || aver
    ngrid = sparse(binx, 1, ones(N, 1), Nf, 1); % no. of obs per grid cell
    if(n0 < 0), n0 = -n0 * N; end % n0 is given as  percentage

    if n0~=0
        % Remove outliers    
        fgrid(ngrid <= n0) = 0;
        ngrid(ngrid <= n0) = 0;
        N = full(sum(ngrid(:))); % how many datapoints are left now?
    end  
end

ind = full(find(fgrid)); % find nonzero values

if aver
    % Make average of values for each point
    fgrid(ind) = fgrid(ind) ./ ngrid(ind);
end

if pad ~= 0
    Nil = find(fgrid == 0);
    fgrid(Nil) = full(fgrid(Nil)) + pad; % Empty grid points are set to pad
end

rho = length(ind) / (Nf); % density of nonzero values in the grid
if rho > 0.4
    fgrid = full(fgrid); 
end
if r == 1
    fgrid = fgrid.';
else
    fgrid = reshape(fgrid, fgsiz);
    switch D % make sure fgrid is stored in the same way as meshgrid
        case 2, fgrid = fgrid.';
        case 3, fgrid = permute(fgrid, [2 1 3]);
    end
end

if (nargout > 0)
    zzgrid = fgrid;
elseif D == 2 % no output, then plot
    colormap(flipud(hot)) %colormap jet
    if 1
        %figure('Position', [100 100 size(fgrid)])    
        imagesc(xvec{:}, fgrid)
        set(gca, 'YDir', 'normal')
    else
        pcolor(xvec{:}, fgrid)
        shading flat %interp
    end
    colorbar
    xlabel(inputname(1))
    ylabel(inputname(1))
    zstr = inputname(2);
    dum = full(size(fgrid'));
    if (~isempty(zstr)) % all this vital information ...
        str = sprintf(['Color scale: %s, %d data points, grid: %dx%d, '...
            'density: %4.2f'], inputname(3), N, dum(1), dum(2), rho);
    else
        str = sprintf('%d data points, grid: %dx%d, density: %4.2f', ...
        N, dum(1), dum(2), rho);
    end
    title(str);
end

end