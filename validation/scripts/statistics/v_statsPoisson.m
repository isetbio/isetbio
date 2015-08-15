%%  v_statsPoisson
%
% Test the iePoisson random number generator calls
%
% BW ISETBIO Team, 2015

ieInit

%% Create a small lambda set of samples and plot them

% A small number
nSamp = 1000;
lambda = 1;
v = iePoisson(lambda,nSamp);
fprintf('Mean %.2f and variance %.2f\n',mean(v(:)), var(v));

vcNewGraphWin;
hist(v,20);


%% A slightly larger number
nSamp = 1000;
lambda = 25;
v = iePoisson(lambda,nSamp);
fprintf('Mean %.2f and variance %.2f\n',mean(v(:)), var(v));

vcNewGraphWin;
hist(v,20);

%%  Now try a big lambda, which should be gaussian

lambda = 1000;
v = iePoisson(lambda,nSamp);
fprintf('Mean %.2f and variance %.2f\n',mean(v(:)), var(v));

vcNewGraphWin;
hist(v,20);


%% END
