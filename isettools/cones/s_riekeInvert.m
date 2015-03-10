%% s_riekeInvert
%    This script computes the physiology differential equation based cone
%    adapted current and use the current to retrieve the input cone
%    absorptions
%
%    Then, we generate some additive Gaussian noise to the adapted current
%    and backward estimate the input cone absorption again. Plots are made
%    to illustrate the effect of cone adaptation noise
%
%  (HJ) ISETBIO TEAM, 2014

%% NOTE: Something is wrong in the program, HJ will fix it soon

%% Compute cone absoprtions and adapted current
% init parameters
s_initISET;
fov = 1;    % field of view
dt = 0.001; % sampling time
expTime = 0.05; % exposure time

% create scene - macbeth color checker
scene = sceneCreate;
scene = sceneSet(scene, 'h fov', fov);
% vcAddObject(scene); sceneWindow;

% create human optics
oi = oiCreate('wvf human');
oi = oiCompute(scene, oi);
% vcAddObject(oi); oiWindow;

% create human sensor
sensor = sensorCreate('human');
sensor = sensorSetSizeToFOV(sensor, fov, scene, oi);
sensor = sensorSet(sensor, 'exp time', dt);

% set eye movement - no eye movement
sensor = sensorSet(sensor, 'sample time interval', dt);
sensor = sensorSet(sensor, 'sensorpositions', zeros(1000, 2));

% compute cone absorptions
sensor = coneAbsorptions(sensor, oi);

% integrate over time
sensor = sensorSet(sensor, 'exp time', expTime);
ph = sensorGet(sensor, 'photons'); ph = cumsum(ph, 3);
ph = ph(:, :, expTime/dt + 1:end) - ph(:, :, 1:end-expTime/dt);
sensor = sensorSet(sensor, 'photons', ph);

% comptue adapted current
[~, adaptedCur] = coneAdapt(sensor, 'rieke');

%% Estimate cone absorption from noise-free adapted current
%  Init parameters for cone adaptation
adaptedCur = - adaptedCur;
p = riekeInit; bgCur = adaptedCur(:,:,1);

% Compute initial state isomerization rate
R = zeros(size(adaptedCur));
bgR = p.sigma * (p.phi*p.smax./(1+(p.q*bgCur/p.beta/p.kGc).^p.n) ...
    .* (p.k*p.beta*p.cdark./(p.beta*p.cdark+p.q*bgCur)./bgCur).^(1/p.h) ...
        - p.eta);
R(:,:,1) = bgR;

% Compute status of intermediate chemicals at initial state
opsin = bgR / p.sigma; PDE = (opsin + p.eta) / p.phi;
Ca  = bgCur * p.q / p.beta; Ca_slow = Ca;
st = p.smax ./ (1 + (Ca / p.kGc).^p.n);
cGMP = st * p.phi ./ (opsin + p.eta);

% Estimate cone absorptions from adapted current
for t = 2 : size(adaptedCur, 3)
    cGMP_old = cGMP; PDE_old = PDE; opsin_old = opsin;
    Ca = Ca + dt * (p.q * adaptedCur(:,:,t-1) - p.beta * Ca);
    Ca_slow = Ca_slow + dt * p.betaSlow * (Ca_slow - Ca);
    cGMP = (adaptedCur(:,:,t) .* (1 + Ca_slow / p.cdark)/p.k).^(1/p.h);
    st = p.smax ./ (1 + (Ca/p.kGc).^p.n);
    PDE = (st - (cGMP - cGMP_old)/dt) ./ cGMP_old;
    opsin = (PDE - PDE_old)/dt  + p.phi * PDE_old - p.eta;
    R(:,:,t) = (opsin - opsin_old)/dt + p.sigma * opsin_old;
end

% Compare estimated absorptions and its true value
indx = (ph > 20); errR = abs(R * expTime - ph) ./ ph;
fprintf('Average error rate:%.2f%%\n', 100 * mean(errR(indx)));

%% Estimate cone absorption from noisy adapted current