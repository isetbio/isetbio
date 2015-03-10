%% v_coneAdaptation
%
%    This script makes comparison between steady state cone adapt
%    computation (see coneAdapt, case 3) and physiological differential
%    equations
%
%  (HJ) ISETBIO TEAM, 2014

%% Init parameters for physiological differential equation
%  These parameters are copied from Fred's script. Sadly, HJ cannot find
%  the original value in the literatures

sigma = 100;        % rhodopsin activity decay rate (1/sec)
phi = 50;			% phosphodiesterase activity decay rate (1/sec)
eta = 100;			% phosphodiesterase activation rate constant (1/sec)
gdark = 35;			% concentration of cGMP in darkness
k = 0.02;       	% constant relating cGMP to current
cdark = 0.5;		% dark calcium concentration
beta = 50;			% rate constant for calcium removal in 1/sec
betaSlow = 2;
n = 4;  			% cooperativity, hill coef
kGc = 0.35;         % hill affinity
h = 3;

q    = 2 * beta * cdark / (k * gdark^h);  % get q using steady state
smax = eta/phi * gdark * (1 + (cdark / kGc)^n);

%% Comparison of steady state response
%  The stimulus used here is static stimulus. We compare the normalized
%  cone absorption response computed by physiological differential
%  equations and Felice A. Dunn et al. model for isomerization rate between
%  1 and 1e5 R*/sec.

%  Compute steady state response for Felice A. Dunn model
rRange  = logspace(0, 5, 2000);
rFelice = 1./ (1 + (45000./ rRange).^0.7);

%  Compute steady state response for physiological differential equations
rPDE = zeros(length(rRange), 1);
for ii = 1 : length(rRange)
    r = rRange(ii);
    fx = @(x) abs(x - k*beta*cdark*smax^h * phi^h /(r/sigma + eta)^h ...
                / (beta*cdark + q*x)/(1 + (q*x/beta/kGc)^n)^h);
    rPDE(ii) = fminbnd(fx, 0, 1000);
end

%  Plot
figure; hold on; grid on;
plot(rRange, rFelice, 'r', 'lineWidth', 2);
plot(rRange, 1./(1+45000./rRange), 'g', 'lineWidth', 2);
plot(rRange, (1 - rPDE/max(rPDE)), 'b', 'lineWidth', 2);
xlabel('isomerization rate (P*/sec)');
ylabel('Nomalized adapted response');
set(gca, 'xscale', 'log');
legend('Dunn Model', 'Linear reciprocal', 'Differential Equation');
title('Steady state response for different models (cone adapt)');

%% Peak ratio of dark adapted and steady state stimulus
%  This computes the ratio of peaks of 1) step stimulus 2) steady state
%  stimulus
%  This part tries to replicate figure 1b (or 2f) for paper:
%       Felice A. Dunn, et al. Light adaptation in cone vision involves
%       switching between receptor and post-receptor sites, 
%       doi:10.1038/nature06150

%  Compute ratio by Dunn's model
ampRatioFelice = (100 + 1.3 * rRange)./(100+rRange)./(0.00029*rRange + 1);
ampRatioFelice(ampRatioFelice > 1) = 1;

% Compute ratio by psysiological differential equations
% init parameters
dt      = 0.001;  % discretize time step - 1 ms
numPts  = 2000;   % total simulation time

opsin   = repmat(rRange/sigma, [numPts 1]);
PDE     = (opsin + eta) / phi;
Ca      = repmat(rPDE'*q/beta, [numPts, 1]);
Ca_slow = Ca;
st      = smax ./ (1 + (Ca / kGc).^n);
cGMP    = st * phi ./ (opsin + eta);

% simulate differential equations
stimulus = 23000;
for ii = 2 : numPts
    opsin(ii, :) = opsin(ii-1, :) + dt*(stimulus - sigma * opsin(ii-1,:));
	PDE(ii, :) = PDE(ii-1, :) + dt * (opsin(ii-1, :)+eta-phi*PDE(ii-1, :));
	Ca(ii, :)  = Ca(ii-1, :) + dt * (q * k * cGMP(ii-1, :).^h ./ ...
                    (1 + Ca_slow(ii-1, :)/cdark) - beta * Ca(ii-1, :));
	Ca_slow(ii, :) = Ca_slow(ii-1, :) - ...
                        dt * betaSlow * (Ca_slow(ii-1, :) - Ca(ii-1, :));
	st(ii, :) = smax ./ (1 + (Ca(ii, :) / kGc).^n);
	cGMP(ii, :) = cGMP(ii-1, :) + ...
                        dt * (st(ii-1, :) - PDE(ii-1, :) .* cGMP(ii-1, :));
end

% compute current
cur = k * cGMP.^h ./ (1 + Ca_slow / cdark);

% compute ampRatio
ampRatioPDE = (rPDE' - min(cur)) ./ (beta*cdark/q - min(cur(:,1)));

% Plot
figure; hold on; grid on;
plot(rRange, ampRatioFelice, 'r', 'lineWidth', 2);
plot(rRange, ampRatioPDE, 'b', 'lineWidth', 2);
set(gca, 'xScale', 'log');
xlim([1 stimulus]); ylim([0 1]);
xlabel('background intensity (P*/sec)');
ylabel('Amp/Amp_{dark}');
legend('Dunn Model', 'Differential Equation');
title('Background dependence of response amplitude');