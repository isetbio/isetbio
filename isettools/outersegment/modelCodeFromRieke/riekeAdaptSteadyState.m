function p = riekeAdaptSteadyState(bgR, p, varargin)
% Steady-state background current calculated from the background rates.
%
%    initialState = riekeAdaptSteadyState(bgR,p,sz)
%
% Inputs:
%  bgR:  Vector (or matrix) of background isomerization rates
%  p:    Parameter list (from riekeInit)
%
%  varargin{1}: sz - sensor array size, (e.g., sensorGet(sensor,'size')),
%                    this parameter is only useful when input bgR is a
%                    scalar
%
% Output:
%  initialState:  The parameters in p augmented by additional terms needed
%                 for the dynamic calculation in riekeAdaptTemporal
%
% Example:
%   riekeAdaptSteadyState(1000,[],[1 1])
%
% The derivation for the formula can be found from:
%   
%
% HJ/VISTASOFT Team, 2014

%% Programming note
%
% Notice that the computation is a search over a bounded variable.  The
% upper and lower bounds are huge for current and thus good enough to
% always find something, we think.

%% Parameters
if notDefined('bgR'), error('Background isomerization rate required.'); end
if isscalar(bgR) && ~isempty(varargin), bgR = bgR*ones(varargin{1});  end
if notDefined('p'),  p = riekeInit; end

%% Calculation
%  In most cases, the input bgR matrix will only contain only a small
%  number of unique values (the background is almost flat) and we can just
%  loop over them
sz = size(bgR);
[bgR, ~, recover_index] = unique(bgR(:));

% Init steady state current for background
bgCur = zeros(length(bgR));
for ii = 1 : length(bgR)
    v = bgR(ii);
    bgCur(ii) = fminbnd(@(x) abs(x - (p.k*p.beta*p.cdark) * ...
        (p.smax*p.phi)^p.h / (v/p.sigma + p.eta)^p.h / ...
        (p.beta*p.cdark + p.q*x) / (1 + (p.q*x/p.beta/p.kGc)^p.n)^p.h), ...
        0, 1000);
end

% recover background isomerization rate and adapted steady state current to
% the input matrix format
bgR   = reshape(bgR(recover_index), sz);
bgCur = reshape(bgCur(recover_index), sz);

% Compute additional initial values
p.opsin   = bgR / p.sigma;
p.PDE     = (p.opsin + p.eta) / p.phi;
p.Ca      = bgCur * p.q / p.beta;
p.Ca_slow = p.Ca;
p.st      = p.smax ./ (1 + (p.Ca / p.kGc).^p.n);
p.cGMP    = p.st * p.phi ./ (p.opsin + p.eta);

p.bgCur = bgCur;

end