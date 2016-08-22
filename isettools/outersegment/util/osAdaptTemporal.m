 function [adaptedData, obj] = osAdaptTemporal(pRate,obj)
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
% if ~exist('p','var') || isempty(p), error('Initial state required.'); end
    
% One millisecond time step
if isfield(obj, 'timeInterval'), dt = obj.timeInterval; else dt = 0.001; end

%% Simulate differential equations
adaptedData = zeros([size(obj.opsin) size(pRate, 3)+1]);
adaptedData(:,:,1) = obj.bgCur;
for ii = 1 : size(pRate, 3)
    obj.opsin = obj.opsin + dt*(obj.OpsinGain*pRate(:,:,ii) - obj.sigma*obj.opsin);
    obj.PDE   = obj.PDE + dt*(obj.opsin + obj.eta - obj.phi * obj.PDE);
    obj.Ca    = obj.Ca + dt*(obj.q*obj.k * obj.cGMP.^obj.h./(1+obj.Ca_slow/obj.cdark)-obj.beta*obj.Ca);
    obj.Ca_slow = obj.Ca_slow - dt * obj.betaSlow * (obj.Ca_slow - obj.Ca);
    obj.st    = obj.smax ./ (1 + (obj.Ca / obj.kGc).^obj.n);
    obj.cGMP  = obj.cGMP  + dt * (obj.st - obj.PDE .* obj.cGMP);
    
    adaptedData(:,:,ii) = - obj.k * obj.cGMP.^obj.h ./ (1 + obj.Ca_slow / obj.cdark);
end

adaptedData(:, :, size(pRate, 3)+1) = adaptedData(:, :, size(pRate, 3));
adaptedData = adaptedData(:,:, 2:end);

end
