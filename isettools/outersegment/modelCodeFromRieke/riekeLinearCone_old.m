function adaptedData = riekeLinearCone_old(pRate,p)

%%
if ~exist('pRate','var') || isempty(pRate), error('Photon absorption rate required.'); end
if ~exist('p','var') || isempty(p), error('Initial state required.'); end
    
% One millisecond time step, as in s_fredConeModel.m
if isfield(p, 'timeInterval'), dt = p.timeInterval; else dt = 0.001; end

%% cone IRF

TimeAxis = [1:size(pRate, 3)]*dt;

ScFact = 0.6745; % To get amplitude right
TauR = 0.0216;   % Rising Phase Time Constant
TauD = 0.0299;   % Damping Time Constant
TauP = 0.5311;   % Period

Io = 2250;                     % half-desensitizing background (in R*/cone/sec, from Juan's paper - corrected)
Ib = pRate(1);
gain_dark = 0.32;              % from Juan's paper (approximate peak of the IRF measured in darkness, and in units of pA/R*) - corrected
gainRatio = 1 ./ (1+(Ib./Io)); % the right side of the equation above, and the gain ratio implied by the bkgnd adapting field

Phi = 34.1814;   % Phase

Filter = ScFact .* (((TimeAxis./TauR).^3)./(1+((TimeAxis./TauR).^3))) .* exp(-((TimeAxis./TauD))).*cos(((2.*pi.*TimeAxis)./TauP)+(2*pi*Phi/360));

newGain = gainRatio .* gain_dark;
oldGain = max(Filter);
Filter = Filter * dt * newGain / oldGain;

for x = 1:size(pRate, 1)
    for y = 1:size(pRate, 2)
        tempData = conv(squeeze(pRate(x, y, :)), Filter);
        adaptedData(x, y, :) = tempData(1:length(Filter)); 
    end
end

end
