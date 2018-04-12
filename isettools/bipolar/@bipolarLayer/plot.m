function hdl = plot(obj, pType, varargin)
% Plot the values from one of the mosaics of a bipolarLayer object
%
% Syntax:
%   hdl = obj.plot(pType, ...)
%
% Description:
%    An interface to layer and mosaic plotting.  At some point, we
%    will have plots that span multiple layers.  For now, this routine
%    mostly reads the mosaic number and passes the request on to
%    that mosaic.
%
% Inputs:
%    pType     - The plot type. Use @bipolarLayer.plot('help','nMosaic',1)
%                to see implemented plot types.
%
% Outputs:
%    hdl       - bipolar mosaic plot according to the provided parameters
%
% Optional Key/Value Pairs:
%    nMosaic   - Mosaic number to plot   (0.  Throws an error)
%    gamma     - controls image display  (1)
%    pos       - positions to plot for time series [1,1]
%    newWindow - whether or not to open in a new window (false)
%
% Notes:

% History:
%    xx/xx/xx jrg/bw (c) isetbio team
%    10/17/17  jnm  Comments & Formatting
%    12/21/17  BW   Cleared the notes by adding defaults and other
%                   comments.

% Examples:
%{
   % ETTBSkip - Example is broken. Variables undefined. Remove this line
   % when fixed.
   s_initRetina;
   bpL.plot('help')
   bpL.plot('mosaic', 'nMosaic', 1);
   bpL.plot('response movie', 'nMosaic', 1);
%}

%% Parse inputs
p = inputParser;
p.CaseSensitive = false;
p.FunctionName  = mfilename;
p.KeepUnmatched = true;   % Sent on to mosaic plot

p.addRequired('pType', @ischar);

p.addParameter('gamma', 1, @isscalar);
p.addParameter('pos', [1, 1], @isvector);
p.addParameter('newWindow', false, @islogical);

% Will this be one of the mosaics, or use multiple mosaics (0, or maybe a
% vector?)
p.addParameter('nMosaic', 0, @(x)(isscalar(x) && (x >= 0)));

% Parse pType
p.parse(pType, varargin{:});

nMosaic = p.Results.nMosaic;
if nMosaic > length(obj.mosaic)
    error('nMosaic (%d) exceeds number of bipolar layers (%d)\n', ...
        nMosaic, length(obj.mosaic));
end

%% Account for parameters
if p.Results.newWindow; hdl = vcNewGraphWin; end

if nMosaic == 0
    % A plot that uses more than one mosaic. What the layer is for.
    % We need to make some stuff up for here.
    disp('No layer plots implemented. This passes to mosaic plot.')
    return;
else
    % A plot based on one mosaic. Call the bipolarMosaic.plot funciton.
    hdl = obj.mosaic{nMosaic}.plot(pType, p.Results);
end
end
