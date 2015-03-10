function p = riekeInit
% Initialize parameters in Rieke adaptation model
%
%    p = riekeInit
%
% See also:  riekeAdaptSteadyState, riekeAdaptTemporal, coneAdapt
%  
% HJ, ISETBIO Team Copyright 2014

% Init parameters
sigma = 100;  % rhodopsin activity decay rate (1/sec)
phi = 50;     % phosphodiesterase activity decay rate (1/sec)
eta = 100;	  % phosphodiesterase activation rate constant (1/sec)
gdark = 35;	  % concentration of cGMP in darkness
k = 0.01;     % constant relating cGMP to current
cdark = 0.5;  % dark calcium concentration
beta = 50;	  % rate constant for calcium removal in 1/sec
betaSlow = 2;
n = 4;  	  % cooperativity, hill coef
kGc = 0.35;   % hill affinity
h = 3;

% Compute more parameters
q    = 2 * beta * cdark / (k * gdark^h);
smax = eta/phi * gdark * (1 + (cdark / kGc)^n);

% Used for finding the steady-state current
p = struct('sigma',sigma, 'phi',phi, 'eta',eta, 'gdark',gdark,...
    'k',k,'cdark',cdark,'beta',beta,'betaSlow',betaSlow,  ...
    'n',n,'kGc',kGc,'h',h,'q',q,'smax',smax');

end