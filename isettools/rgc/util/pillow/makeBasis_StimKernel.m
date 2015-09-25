function [kbas,kbasis] = makeBasis_StimKernel(kbasprs, nkt);
%  [kbas, kbasis] = makeBasis_StimKernel(kbasprs, nkt);
% 
%  Generates a basis consisting of raised cosines and several columns of
%  identity matrix vectors for temporal structure of stimulus kernel
%
%  Args: kbasprs = struct with fields: 
%          neye = number of identity basis vectors at front
%          ncos = # of vectors that are raised cosines
%          kpeaks = 2-vector, with peak position of 1st and last vector,
%             relative to start of cosine basis vectors (e.g. [0 10])
%          b = offset for nonlinear scaling.  larger values -> more linear
%             scaling of vectors.  bk must be >= 0
%        nkt = number of time samples in basis (optional)
%
%  Output:
%        kbas = orthogonal basis
%        kbasis = standard (non-orth) basis

neye = kbasprs.neye;
ncos = kbasprs.ncos;
kpeaks = kbasprs.kpeaks;
b = kbasprs.b;


kdt = 1;  % spacing of x axis must be in units of 1

% nonlinearity for stretching x axis (and its inverse)
nlin = @(x)log(x+1e-20);
invnl = @(x)exp(x)-1e-20; % inverse nonlinearity


% Generate basis of raised cosines
yrnge = nlin(kpeaks+b);  
db = diff(yrnge)/(ncos-1);      % spacing between raised cosine peaks
ctrs = yrnge(1):db:yrnge(2);  % centers for basis vectors
mxt = invnl(yrnge(2)+2*db)-b; % maximum time bin
kt0 = [0:kdt:mxt]';
nt = length(kt0);        % number of points in iht
ff = @(x,c,dc)(cos(max(-pi,min(pi,(x-c)*pi/dc/2)))+1)/2; % raised cosine basis vector
kbasis0 = ff(repmat(nlin(kt0+b), 1, ncos), repmat(ctrs, nt, 1), db);

% Concatenate identity-vectors
nkt0 = size(kt0,1);
kbasis = [[eye(neye); zeros(nkt0,neye)] [zeros(neye, ncos); kbasis0]];
kbasis = flipud(kbasis);  % flip so fine timescales are at the end.
nkt0 = size(kbasis,1);

if nargin > 1
    if nkt0 < nkt 
        %fprintf('>>>> makeBasis_StimKernel: Padding basis with %d rows of zeros\n', nkt-nkt0);
        kbasis = [zeros(nkt-nkt0,ncos+neye); kbasis];
    elseif nkt0 > nkt
        %fprintf('>>>> makeBasis_StimKernel: Removing %d rows from basis\n', nkt0-nkt);
        kbasis = kbasis(end-nkt+1:end,:);
    end
end

kbasis = normalizecols(kbasis);
%kbas = orth(kbasis);
kbas = kbasis;
