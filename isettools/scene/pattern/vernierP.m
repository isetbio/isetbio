function params = vernierP(varargin)
% Default vernier stimulus parameters
%
%   p = vernierP;
%
% To specify values for the params use param/val pairs
%
% BW, ISETBIO Team, 2016

%% Parse arguments

p = inputParser;

p.addParameter('display',displayCreate('LCD-Apple'),@isstruct);
p.addParameter('sceneSz',[50 50],@isnumeric);
p.addParameter('name','unknown',@ischar);

p.addParameter('offset',1,@isinteger);
p.addParameter('bgColor',0.5,@isscalar);

p.addParameter('barWidth',1,@isnumeric);
p.addParameter('barLength',[],@isscalar);
p.addParameter('barColor',1,@isscalar);

p.parse(varargin{:});

%% Assign parameters

% Identifier
params.name  = p.Results.name;           % Char

% Display scene
params.display  = p.Results.display;     % Display structure
params.sceneSz  = p.Results.sceneSz;     % Pixels

% General
params.offset   = p.Results.offset;      % Pixels
params.bgColor  = p.Results.bgColor;     % 0 to 1?

% Bar properties.  Define better
params.barWidth  = p.Results.barWidth;    % Pixels
params.barLength = p.Results.barLength;   % Pixels?
params.barColor  = p.Results.barColor;    % (rgb?)

params.name      = p.Results.barColor;    % Identifier

end
