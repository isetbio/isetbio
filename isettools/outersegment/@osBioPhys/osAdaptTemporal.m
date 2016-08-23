 function [adaptedData, model] = osAdaptTemporal(pRate,obj)
% Time varying voltage response from photon rate and initial state
%
%    adaptedData = osAdaptTemporal(pRate, obj)
%                    [only called internally from @osBioPhys/osCompute.m]
% 
% Inputs:
%   pRate - Photon absorption rate from coneMosaic.absorptions/coneMosaic.integrationTime
%   obj   - osBioPhys object containing many initial parameters
%
% Return
%   adaptedData - adapted photocurrent data (pA) for coneMosaic.current
%   obj         - osBioPhys object containing many final parameters
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
% This model gives a cone-by-cone adaptation and produces a time-series
% structure in adaptedDat that is stored into the current field of the cone
% mosaic object in @osBioPhys/osCompute.m.
% 
% Reference:
%   http://isetbio.org/cones/adaptation%20model%20-%20rieke.pdf
%   https://github.com/isetbio/isetbio/wiki/Cone-Adaptation
%
% Example:
%   From @osBioPhys/osCompute.m, line 64:
%        [current, model.state]  = osAdaptTemporal(pRate, model.state);
%
% See also:
%   osAdaptSteadyState, osAdaptTemporal
%
% 
% HJ, ISETBIO Team, 2014
% JRG, ISETBIO Team, updated 8/2016

%%  Check inputs
if ~exist('pRate','var') || isempty(pRate), error('Photon absorption rate required.'); end
    
dt = obj.timeStep;
model = obj.state;
%% Simulate differential equations
adaptedData = zeros([size(model.opsin) size(pRate, 3)+1]);
adaptedData(:,:,1) = model.bgCur;

q    = 2 * model.beta * model.cdark / (model.k * model.gdark^model.h);
smax = model.eta/model.phi * model.gdark * (1 + (model.cdark / model.kGc)^model.n);  

for ii = 1 : size(pRate, 3)
    model.opsin = model.opsin + dt*(model.OpsinGain*pRate(:,:,ii) - model.sigma*model.opsin);
    model.PDE   = model.PDE + dt*(model.opsin + model.eta - model.phi * model.PDE);
    model.Ca    = model.Ca + dt*(q*model.k * model.cGMP.^model.h./(1+model.Ca_slow/model.cdark)-model.beta*model.Ca);
    model.Ca_slow = model.Ca_slow - dt * model.betaSlow * (model.Ca_slow - model.Ca);
    model.st    = smax ./ (1 + (model.Ca / model.kGc).^model.n);
    model.cGMP  = model.cGMP  + dt * (model.st - model.PDE .* model.cGMP);
    
    adaptedData(:,:,ii) = - model.k * model.cGMP.^model.h ./ (1 + model.Ca_slow / model.cdark);
end

adaptedData(:, :, size(pRate, 3)+1) = adaptedData(:, :, size(pRate, 3));
adaptedData = adaptedData(:,:, 2:end);
model = obj.model;
end
