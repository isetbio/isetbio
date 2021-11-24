function [ih, ihbasprs] = buildPostSpikeFilter(dt, varargin)
% Builds the post spike filters for the rgcMosaicGLM object.
%
% Syntax:
%   [ih, ihbasprs] = buildPostSpikeFilter(dt, [varargin])
%
% Description:
%    This is the original Pillow code used in the GLM.  We tried not to
%    modify it, and indeed it is now set just as in the original apart from
%    some editing of the comments and parameters for clarity, and to return
%    the parameters for reproducibility.
%
%    We might consider whether we want to add an input argument so that we
%    can set the parameters for experiments.  In that case, we might decide
%    to rename the routine, referencing the original.
%
%    This function contains examples of usage inline. To access these, type
%    'edit buildPostSpikeFilter.m' into the Command Window.
%
% Inputs:
%    dt      - Numeric. The bin fraction of astimulus time step. The value
%              for dt is typically <= 0.01. What are the units? sec? ms?
%
% Outputs: 
%    ih      - Array. The postSpikeFilters as a cell array, where each cell
%              is a time series filters. Note: these are not "filters" in
%              the convolutional sense; in computeSpikesGLM, they are added
%              to the time series response of the particular neuron.
%
% Optional key/value pairs:
%    nCols   - Numeric. The number of basis vectors for the post-spike
%              kernel. Default is 5.
%    hpeaks  - Array. An array of length two containing the peak locations
%              for the first and last vectors. Default [.1 2].
%    b       - Numeric. The measurement for how nonlinear to make the
%              spacings. Default 0.5.
%    absref  - Numeric. The absolute refractory period. Default 0.1.
%    weights - Array. An array of length 5 containing parameter weights.
%              Default [-10 -5 0 2 -2].
%
% Notes:
%    * TODO: Ask EJ why the values in ih are fixed.
%

% History:
%    XX/XX/XX  JP   Written by J. Pillow
%    XX/XX/XX  BW   Wandell - modified for clarity & return the parameters
%    06/10/19  JNM  Documentation pass, completed programming TODO for
%                   ability to override default parameters.

% Examples:
%{
    obj.postSpikeFilter = buildPostSpikeFilter(.01);
%}

% % --- Make basis for post-spike (h) current ------
% ihbasprs.ncols = 5;       % Number of basis vectors for post-spike kernel
% ihbasprs.hpeaks = [.1 2]; % Peak location for first and last vectors
% ihbasprs.b = .5;          % How nonlinear to make spacings
% ihbasprs.absref = .1;     % absolute refractory period 
% ihbasprs.weights = [-10 -5 0 2 -2]';

% [Note: JNM attempting programming TODO.]
p = inputParser;
p.addParameter('nCols', 5, @isscalar);
p.addParameter('hpeaks', [.1 2], ...
    @(x) isnumeric(x) && (isequal(length(x), 2)));
p.addParameter('b', 0.5, @isscalar);
p.addParameter('absref', 0.1, @isscalar);
p.addParameter('weights', [-10 -5 0 2 -2]', ...
    @(x) isnumeric(x) && (isequal(length(x), 5)));
p.parse(varargin{:});

ihbasprs.ncols = p.Results.nCols;
ihbasprs.hpeaks = p.Results.hpeaks;
ihbasprs.b = p.Results.b;
ihbasprs.absref = p.Results.absref;
ihbasprs.weights = p.Results.weights;

% Original
% [iht, ihbas, ihbasis] = makeBasis_PostSpike(ihbasprs, dt);
[~, ~, ihbasis] = makeBasis_PostSpike(ihbasprs, dt);

% Why are these values fixed? Ask EJ.
ih = ihbasis * ihbasprs.weights;  % h current

% Original code
% ph = 1;

end
