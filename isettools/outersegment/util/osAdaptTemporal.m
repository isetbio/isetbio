 function [adaptedData, p] = osAdaptTemporal(pRate,p)
% Time varying voltage response from photon rate and initial state
%
%    adaptedData = osAdaptTemporal(pRate,initialState)
% 
% Inputs:
%   pRate - Photon absorption rate
%   p     - Structure containing many initial parameters
%
% Return
%   adaptedData - adapted photocurrent data (pA)
%   p           - structure of final state parameters
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
%   t_os, t_osIntroduction
%
% See also:
%   osInit, osAdaptSteadyState, coneAdapt
%
% HJ, ISETBIO Team, 2014

%%  Check inputs
if ~exist('pRate','var') || isempty(pRate), error('Photon absorption rate required.'); end
if ~exist('p','var') || isempty(p), error('Initial state required.'); end
    
% One millisecond time step
if isfield(p, 'timeInterval'), dt = p.timeInterval; else dt = 0.001; end

%% Simulate differential equations
adaptedData = zeros([size(p.opsin) size(pRate, 3)+1]);
adaptedData(:,:,1) = p.bgCur;
for ii = 1 : size(pRate, 3)
    p.opsin = p.opsin + dt*(p.OpsinGain*pRate(:,:,ii) - p.sigma*p.opsin);
    p.PDE   = p.PDE + dt*(p.opsin + p.eta - p.phi * p.PDE);
    p.Ca    = p.Ca + dt*(p.q*p.k * p.cGMP.^p.h./(1+p.Ca_slow/p.cdark)-p.beta*p.Ca);
    p.Ca_slow = p.Ca_slow - dt * p.betaSlow * (p.Ca_slow - p.Ca);
    p.st    = p.smax ./ (1 + (p.Ca / p.kGc).^p.n);
    p.cGMP  = p.cGMP  + dt * (p.st - p.PDE .* p.cGMP);
    
    adaptedData(:,:,ii) = - p.k * p.cGMP.^p.h ./ (1 + p.Ca_slow / p.cdark);
end

adaptedData(:, :, size(pRate, 3)+1) = adaptedData(:, :, size(pRate, 3));
adaptedData = adaptedData(:,:, 2:end);

end
