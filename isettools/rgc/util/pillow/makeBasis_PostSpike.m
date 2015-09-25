function [iht, ihbas, ihbasis] = makeBasis_PostSpike(ihprs,dt);
% [iht, ihbas, ihbasis] = makeBasis_PostSpike(ihprs,dt);
%
% Make nonlinearly stretched basis consisting of raised cosines
% Inputs: prs: param structure with fields:
%          ncols = # of basis vectors
%          hpeaks = 2-vector containg [1st_peak  last_peak], the peak 
%                      location of first and last raised cosine basis vectors
%          b = offset for nonlinear stretching of x axis:  y = log(x+b) 
%                 (larger b -> more nearly linear stretching)
%          absref = absolute refractory period (optional)
%
%  Outputs:  iht = time lattice on which basis is defined
%            ihbas = orthogonalized basis
%            ihbasis = original (non-orthogonal) basis 
%
%  Example call:
%
%  ihbasprs.ncols = 5;  
%  ihbasprs.hpeaks = [.1 2];  
%  ihbasprs.b = .5;  
%  ihbasprs.absref = .1;  %% (optional)
%  [iht,ihbas,ihbasis] = makeBasis_PostSpike(ihprs,dt);

ncols = ihprs.ncols;
b = ihprs.b;
hpeaks = ihprs.hpeaks;
if isfield(ihprs, 'absref');
    absref = ihprs.absref;
else
    absref = 0;
end

% Check input values
if (hpeaks(1)+b) < 0, 
    error('b + first peak location: must be greater than 0'); 
end
if absref >= dt  % use one fewer "cosine-shaped basis vector
    ncols = ncols-1;
elseif absref > 0
    warning('Refractory period too small for time-bin sizes');
end

% nonlinearity for stretching x axis (and its inverse)
nlin = @(x)log(x+1e-20);
invnl = @(x)exp(x)-1e-20; % inverse nonlinearity

% Generate basis of raised cosines
yrnge = nlin(hpeaks+b);  
db = diff(yrnge)/(ncols-1);    % spacing between raised cosine peaks
ctrs = yrnge(1):db:yrnge(2);   % centers for basis vectors
mxt = invnl(yrnge(2)+2*db)-b;  % maximum time bin
iht = [0:dt:mxt]';
nt = length(iht);        % number of points in iht
ff = @(x,c,dc)(cos(max(-pi,min(pi,(x-c)*pi/dc/2)))+1)/2; % raised cosine basis vector
ihbasis = ff(repmat(nlin(iht+b), 1, ncols), repmat(ctrs, nt, 1), db);

% set first cosine basis vector bins (before 1st peak) to 1
ii = find(iht<=hpeaks(1));
ihbasis(ii,1) = 1;

% create first basis vector as step-function for absolute refractory period
if absref >= dt
    ii = find(iht<absref);
    ih0 = zeros(size(ihbasis,1),1);
    ih0(ii) = 1;
    ihbasis(ii,:) = 0;
    ihbasis = [ih0,ihbasis];
end
ihbas = orth(ihbasis);  % use orthogonalized basis
