function val = lorentzSum(params, x)
%% Compute sum of Lorentzian components
%    val = lorentzSum(params, x)
%
%  This function computes sum of output of multiple lorentzian components.
%  Value of Lorentzian component can be computed as
%    y = S / (1 + (x/f)^2)^n
%  Here, S, f and n are constant parameters
%
%  Inputs:
%    params - n-by-3 parameter matrix, each row contains S,f,n values for
%             one Lorentzian component
%    x      - positions to be evaluated
%
%  Output:
%    val  - sum of output of Lorentzian components, same size as x
%
% (HJ) ISETBIO TEAM, 2014

%% Check inputs
if ~exist('params', 'var'), error('parameters required'); end
if isvector(params), params = params(:)'; end
if size(params, 2) ~= 3, error('parameter matrix size error'); end
if ~exist('x', 'var'), error('evaluation point x required'); end

%% Accumulate values from each component
%  Init val to zero
val = zeros(size(x));

%  Make sure params are non-negative
params = abs(params);

%  Loop and compute val for each component
for ii = 1 : size(params, 1)
    val = val + params(ii, 1) ./ (1+(x/params(ii, 2)).^2).^ params(ii, 3);
end

end