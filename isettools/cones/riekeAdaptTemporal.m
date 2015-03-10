function adaptedData = riekeAdaptTemporal(pRate,p)
% Time varying voltage response from photon rate and initial state
%
%    adaptedData = riekeAdaptTemporal(pRate,initialState)
% 
% Inputs:
%   pRate: Photon absorption rate
%   initialState:  Structure containing many initial parameters
%
% Return
%   adaptedVolts:
%
% In this case, the physiological differential equations for cones are
% implemented. The differential equations are:
%
%    1) d opsin(t) / dt = -sigma * opsin(t) + R*(t)
%    2) d PDE(t) / dt = opsin(t) - phi * PDE(t) + eta
%    3) d cGMP(t) / dt = S(t) - PDE(t) * cGMP(t)
%    4) d Ca(t) / dt = q * I(t) - beta * Ca(t)
%    5) d Ca_slow(t) / dt = - beta_slow * (Ca_slow(t) - Ca(t))
%    6) S(t) = smax / (1 + (Ca(t)/kGc)^n)
%    7) I(t) = k * cGMP(t)^h / (1 + Ca_slow/Ca_dark)
%
% This model gives a cone-by-cone adaptation and it requires a time-series
% data in sensor structure
%
% Example:
%   See coneAdapt???  Or some script somewhere that needs to be written.
%
% See also:  riekeInit, riekeAdaptSteadyState, coneAdapt
%
% HJ ISETBIO Team, Copyright 2014

%%
if ~exist('pRate','var') || isempty(pRate), error('Photon absorption rate required.'); end
if ~exist('p','var') || isempty(p), error('Initial state required.'); end
    
% One millisecond time step, as in s_fredConeModel.m
if isfield(p, 'timeInterval'), dt = p.timeInterval; else dt = 0.001; end

%% Simulate differential equations

adaptedData = zeros([size(p.opsin) size(pRate, 3)+1]);
adaptedData(:,:,1) = p.bgCur;
for ii = 1 : size(pRate, 3)
    p.opsin = p.opsin + dt * (pRate(:,:,ii) - p.sigma * p.opsin);
    p.PDE   = p.PDE   + dt * (p.opsin + p.eta - p.phi * p.PDE);
    p.Ca    = p.Ca    + dt * (p.q*p.k * p.cGMP.^p.h./(1+p.Ca_slow/p.cdark)-p.beta*p.Ca);
    p.Ca_slow = p.Ca_slow - dt * p.betaSlow * (p.Ca_slow - p.Ca);
    p.st    = p.smax ./ (1 + (p.Ca / p.kGc).^p.n);
    p.cGMP  = p.cGMP  + dt * (p.st - p.PDE .* p.cGMP);
    
    adaptedData(:,:,ii) = - p.k * p.cGMP.^p.h ./ (1 + p.Ca_slow / p.cdark);
%     adaptedData(:,:,ii) = - p.k * p.cGMP.^p.h ./ (1 + p.Ca_slow / p.cdark);
%     p.cGMP  = p.cGMP  + dt * (p.st - p.PDE .* p.cGMP);
%     p.st    = p.smax ./ (1 + (p.Ca / p.kGc).^p.n);
%     p.Ca_slow = p.Ca_slow - dt * p.betaSlow * (p.Ca_slow - p.Ca);
%     p.Ca    = p.Ca    + dt * (p.q*adaptedData(:,:,ii-1)-p.beta*p.Ca);
%     p.PDE   = p.PDE   + dt * (p.opsin + p.eta - p.phi * p.PDE);
%     p.opsin = p.opsin + dt * (pRate(:,:,ii) - p.sigma * p.opsin);
end

adaptedData = adaptedData(:,:, 2:end);

end
