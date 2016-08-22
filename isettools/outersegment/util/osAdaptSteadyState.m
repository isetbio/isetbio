function obj = osAdaptSteadyState(obj, bgR, varargin)
% Steady-state background current calculated from the background rates.
%
%    initialState = osAdaptSteadyState(bgR,p,sz)
%
% Inputs:
%  bgR:  Vector (or matrix) of background isomerization rates
%  p:    Parameter list (from osInit)
%
%  varargin{1}: sz - sensor array size, (e.g., sensorGet(sensor,'size')),
%                    this parameter is only useful when input bgR is a
%                    scalar
%
% Output:
%  initialState:  The parameters in p augmented by additional terms needed
%                 for the dynamic calculation in osAdaptTemporal
%
% Example:
%   osAdaptSteadyState(1000,[],[1 1])
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
% if notDefined('p'),  p = osInit; end

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
    bgCur(ii) = fminbnd(@(x) abs(x - (obj.k*obj.beta*obj.cdark) * ...
        (obj.smax*obj.phi)^obj.h / (v/obj.sigma + obj.eta)^obj.h / ...
        (obj.beta*obj.cdark + obj.q*x) / (1 + (obj.q*x/obj.beta/obj.kGc)^obj.n)^obj.h), ...
        0, 1000);
end

% recover background isomerization rate and adapted steady state current to
% the input matrix format
bgR   = reshape(bgR(recover_index), sz);
bgCur = reshape(bgCur(recover_index), sz);

% Compute additional initial values
obj.opsin   = bgR / obj.sigma;
obj.PDE     = (obj.opsin + obj.eta) / obj.phi;
obj.Ca      = bgCur * obj.q / obj.beta;
obj.Ca_slow = obj.Ca;
obj.st      = obj.smax ./ (1 + (obj.Ca / obj.kGc).^obj.n);
obj.cGMP    = obj.st * obj.phi ./ (obj.opsin + obj.eta);

obj.bgCur = bgCur;

end