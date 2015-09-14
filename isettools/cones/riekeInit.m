function p = riekeInit
% Initialize parameters in Rieke adaptation model
%
%    p = riekeInit
%
% See also:  riekeAdaptSteadyState, riekeAdaptTemporal, coneAdapt
%  
% HJ, ISETBIO Team Copyright 2014

% Init parameters
sigma = 22;  % rhodopsin activity decay rate (1/sec) - default 22
phi = 22;     % phosphodiesterase activity decay rate (1/sec) - default 22
eta = 2000;	  % phosphodiesterase activation rate constant (1/sec) - default 2000
gdark = 20.5; % concentration of cGMP in darkness - default 20.5
k = 0.02;     % constant relating cGMP to current - default 0.02
h = 3;       % cooperativity for cGMP->current - default 3
cdark = 1;  % dark calcium concentration - default 1
beta = 9;	  % rate constant for calcium removal in 1/sec - default 9
betaSlow = 0.4; % rate constant for slow calcium modulation of channels - default 0.4
n = 4;  	  % cooperativity for cyclase, hill coef - default 4
kGc = 0.5;   % hill affinity for cyclase - default 0.5
OpsinGain = 10; % so stimulus can be in R*/sec (this is rate of increase in opsin activity per R*/sec) - default 10

% Compute more parameters - steady state constraints among parameters
q    = 2 * beta * cdark / (k * gdark^h);
smax = eta/phi * gdark * (1 + (cdark / kGc)^n);

% Used for finding the steady-state current
p = struct('sigma',sigma, 'phi',phi, 'eta',eta, 'gdark',gdark,...
    'k',k,'cdark',cdark,'beta',beta,'betaSlow',betaSlow,  ...
    'n',n,'kGc',kGc,'h',h,'q',q,'smax',smax','OpsinGain',OpsinGain);

end