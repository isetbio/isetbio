%%v_statsPoisson  Test the iePoisson random number generator calls
%
% Description:
%     Exercise the Poisson random number generator

% BW ISETBIO Team, 2015
% 10/19/17  dhb  This was broken because someone changed the calling args
%                to iePoisson without fixing this.  Fixed so it runs again.

%% Initialize
ieInit

%% Create a small lambda set of samples and plot them

% A small number
nSamp = 1000;
lambda = 1;
v = iePoisson(lambda,'nSamp',nSamp);
% fprintf('Mean %.2f and variance %.2f\n',mean(v(:)), var(v));

vcNewGraphWin;
hist(v,20);

%% A slightly larger number
nSamp = 1000;
lambda = 25;
v = iePoisson(lambda,'nSamp',nSamp);
% fprintf('Mean %.2f and variance %.2f\n',mean(v(:)), var(v));

vcNewGraphWin;
hist(v,20);

%%  Now try a big lambda, which should be gaussian
lambda = 1000;
v = iePoisson(lambda,'nSamp',nSamp);
% fprintf('Mean %.2f and variance %.2f\n',mean(v(:)), var(v));
vcNewGraphWin;
hist(v,20);

