function [hdl, uData] = plot(obj, pType, varargin)
% Plot function for rgcLayer
%
% Syntax:
%   [hf, uData] = rgcLayer.plot(type, [varargin])
%
% Description:
%    When the plots are directed to a single mosaic, this call forwards
%    them to @rgcMosaic.plot.
%
%    This function contains examples of usage inline. To access these, type
%    'edit @rgcLayer/plot.m' into the Command Window.
%
% Inputs:
%    pType   - String. The type of plot. Type @rgcLayer.plot('help') to see
%              allowable plot types.
%
% Outputs:
%    hf      - Handle. The figure handle.
%    uData   - Matrix. The computed data.
%
% Optional key/value pairs:
%    nMosaic - Numeric. The number indicator of which mosaic to use.
%              Default is 1.
%    hf      - Handle. The figure handle or control structure, the meaning
%              of value is:
%        []: Default. Create the plot in a new figure
%        figure handle: plot in figure specified by hf
%        'none: don't plot
%

% History:
%    XX/XX/16  JRG/HJ/BW  ISETBIO TEAM, 2016
%    06/21/19  JNM        Documentation pass

% Examples:
%{
    % ETTBSkip - skipping broken examples
    %
    % This needs code to define rgcLayer before it could
    % possibly work.
    rgcLayer.plot('mosaic', 'nMosaic', 1);
    rgcLayer.plot('psth', 'nMosaic', 2);
    rgcLayer.plot('psth', 'nMosaic', 2, 'cell', [1 1]);
%}

%%
p = inputParser;
p.KeepUnmatched = true;

% Set up vFunc for checking the pType
pType = ieParamFormat(pType);
p.addRequired('pType', @ischar);  % Type of plot

% figure handle can be a graphics object, the string 'none', or empty
vFunc = @(x)( ischar(x) || isgraphics(x) || isempty(x));
p.addParameter('hf', [], vFunc);

p.addParameter('nMosaic', 1, @isinteger);
p.parse(pType, varargin{:});

% Which mosaic.
nMosaic = p.Results.nMosaic;
if nMosaic > length(obj.mosaic)
    error('nMosaic (%d) exceeds number of bipolar layers (%d)\n', ...
        nMosaic, length(obj.mosaic));
end

%% Account for parameters
uData = [];
hf = p.Results.hf;
if ischar(hf) && strcmp(hf, 'none')
    hf = [];
elseif isempty(hf) && ~strcmpi(pType, 'help')
    hf = vcNewGraphWin;
end

if nMosaic == 0
    % A plot that uses more than one mosaic. What the layer is for.
    % We need to make some stuff up for here.
    disp('No layer plot type are implemented yet.')
    return;
else
    % A plot based on one mosaic. Call the bipolarMosaic.plot funciton.
    [hdl, uData] = obj.mosaic{nMosaic}.plot(pType, 'hf', hf);
end

end
