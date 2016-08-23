function model = osAdaptSteadyState(obj, bgR, varargin)
% Steady-state background current calculated from the background rates.
%
%    obj = osAdaptSteadyState(obj, bgR, sz)
%                    [only called internally from @osBioPhys/osCompute.m]
%
% Inputs:
%  obj - osBioPhys object containing pre-set parameters from obj init
%  bgR - Vector (or matrix) of background isomerization rates
% 
%  varargin{1}: sz - sensor array size, (e.g., sensorGet(sensor,'size')),
%                    this parameter is only useful when input bgR is a
%                    scalar
%
% Output:
%   obj   - osBioPhys object containing all initial parameters needed for
%       computation of adapted current over time series input.
%
% Example:
%   From @osBioPhys/osCompute.m, line 60:
%     obj.model.state = osAdaptSteadyState(obj, bgR, [size(pRate, 1) size(pRate, 2)]);
% 
% Reference:
%   http://isetbio.org/cones/adaptation%20model%20-%20rieke.pdf
%   https://github.com/isetbio/isetbio/wiki/Cone-Adaptation
%
% HJ, ISETBIO Team, 2014
% JRG, ISETBIO Team, updated 8/2016
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

q    = 2 * obj.model.beta * obj.model.cdark / (obj.model.k * obj.model.gdark^obj.model.h);
smax = obj.model.eta/obj.model.phi * obj.model.gdark * (1 + (obj.model.cdark / obj.model.kGc)^obj.model.n);

model = obj.model;

% Init steady state current for background
bgCur = zeros(length(bgR));
for ii = 1 : length(bgR)
    v = bgR(ii);
    bgCur(ii) = fminbnd(@(x) abs(x - (obj.model.k*obj.model.beta*obj.model.cdark) * ...
        (smax*obj.model.phi)^obj.model.h / (v/obj.model.sigma + obj.model.eta)^obj.model.h / ...
        (obj.model.beta*obj.model.cdark + q*x) / (1 + (q*x/obj.model.beta/obj.model.kGc)^obj.model.n)^obj.model.h), ...
        0, 1000);
end

% recover background isomerization rate and adapted steady state current to
% the input matrix format
bgR   = reshape(bgR(recover_index), sz);
bgCur = reshape(bgCur(recover_index), sz);

% Compute additional initial values
model.opsin   = bgR / obj.model.sigma;
model.PDE     = (model.opsin + obj.model.eta) / obj.model.phi;
model.Ca      = bgCur * q / obj.model.beta;
model.Ca_slow = model.Ca;
model.st      = smax ./ (1 + (model.Ca / obj.model.kGc).^obj.model.n);
model.cGMP    = model.st * obj.model.phi ./ (model.opsin + obj.model.eta);

model.bgCur = bgCur;

end