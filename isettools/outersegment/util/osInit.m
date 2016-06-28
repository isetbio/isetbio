function p = osInit
% Initialize parameters in Rieke adaptation model
%
%    p = osInit
%
% See also:  
%   osAdaptSteadyState, osAdaptTemporal, coneAdapt
%  
% HJ, ISETBIO Team, 2014

% Init parameters
sigma = 22;     % rhodopsin activity decay rate (1/sec)
phi = 22;       % phosphodiesterase activity decay rate (1/sec)
eta = 2000;	    % phosphodiesterase activation rate constant (1/sec)
gdark = 20.5;   % concentration of cGMP in darkness
k = 0.02;       % constant relating cGMP to current
h = 3;          % cooperativity for cGMP->current
cdark = 1;      % dark calcium concentration
beta = 9;	    % rate constant for calcium removal in 1/sec
betaSlow = 0.4; % rate constant for slow calcium modulation of channels
n = 4;  	    % cooperativity for cyclase, hill coef
kGc = 0.5;      % hill affinity for cyclase
OpsinGain = 10; % rate of increase in opsin activity per R*/sec

% Compute more parameters - steady state constraints among parameters
q    = 2 * beta * cdark / (k * gdark^h);
smax = eta/phi * gdark * (1 + (cdark / kGc)^n);

% Used for finding the steady-state current
p = struct('sigma',sigma, 'phi',phi, 'eta',eta, 'gdark',gdark,...
    'k',k,'cdark',cdark,'beta',beta,'betaSlow',betaSlow,  ...
    'n',n,'kGc',kGc,'h',h,'q',q,'smax',smax','OpsinGain',OpsinGain);

end