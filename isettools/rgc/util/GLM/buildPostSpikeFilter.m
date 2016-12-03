function [ih, ihbasprs] = buildPostSpikeFilter(dt)
% Builds the post spike filters for the rgcMosaicGLM object.
% 
%  [ih, ihbasprs] = buildPostSpikeFilter(dt)
%
% This is the original Pillow code used in the GLM.  We tried not to modify
% it, and indeed it is now set just as in the original apart from some
% editing of the comments and parameters for clarity, and to return the
% parameters for reproducibility.
%
% We might consider whether we want to add an input argument so that we can
% set the parameters for experiments.  In that case, we might decide to
% rename the routine, referencing the original.
%
% Inputs: 
%   dt - for bin fraction of stimulus time step
%   (usually 0.01 or less)
%         What are the units? Sec?  Milliseconds? 
%   We should add the option of sending in ihbaspars as an argument and
%   over-riding the parameters within here by those parameters.
% 
% Outputs: 
%   ih - the postSpikeFilters as a cell array, where each cell is a time
%   series filters. Note: these are not "filters" in the convolutional
%   sense; in computeSpikesGLM, they are added to the time series response
%   of the particular neuron.
% 
% Example:
%     obj.postSpikeFilter = buildPostSpikeFilter(.01);
% 
%%%% Written by J. Pillow
%
%
% Modified for clarity and to return the parameters
% B. Wandell

% --- Make basis for post-spike (h) current ------
ihbasprs.ncols = 5;        % Number of basis vectors for post-spike kernel
ihbasprs.hpeaks = [.1 2];  % Peak location for first and last vectors
ihbasprs.b = .5;           % How nonlinear to make spacings
ihbasprs.absref = .1;      % absolute refractory period 
ihbasprs.weights = [-10 -5 0 2 -2]';

% Original
% [iht,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs,dt);
[~,~,ihbasis] = makeBasis_PostSpike(ihbasprs,dt);

% Why are these values fixed?  Ask EJ.
ih = ihbasis*ihbasprs.weights;  % h current

% Original code
% ph=1;

end

