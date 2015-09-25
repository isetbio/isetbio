function ih = buildPostSpikeFilter(dt)

%%%% Written by J. Pillow

% dt = .01
% --- Make basis for post-spike (h) current ------
ihbasprs.ncols = 5;  % Number of basis vectors for post-spike kernel
ihbasprs.hpeaks = [.1 2];  % Peak location for first and last vectors
ihbasprs.b = .5;  % How nonlinear to make spacings
ihbasprs.absref = .1; % absolute refractory period 
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs,dt);
ih = ihbasis*[-10 -5 0 2 -2]';  % h current
ph=1;
