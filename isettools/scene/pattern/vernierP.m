function params = vernierP(varargin)
% Default vernier stimulus parameters
%
% Syntax:
%	p = vernierP([varargin]);
%
% Description:
%    Parameters to specify a vernier scene.  
%
%    Generally  used with 
%         scene = sceneCreate('vernier', 'display', params),
%      or
%         oisCreate('vernier', ....)
%
%    There are examples in the code. To access, type 'edit vernierP.m' into
%    the Command Window.
%
% Inputs:
%    varargin - (Optional) A structure containing the additional parameters
%               to instill in the vernier stimulus. Default is to not
%               provide any, and use the defaults listed with the
%               parameters below. Parameters such as:
%                    name      - Used by oiSequence.visualize.
%                                Default 'unknown'
%                    display   - Display structure
%                                Default type 'LCD-Apple' with default wave
%                    sceneSz   - Pixels on display
%                                Default [50 50]
%                    offset    - Pixels on display
%                                Default 1
%                    bgColor   - Background color, e.g., [.6 .4 .2]
%                                Default 0.5
%                    barWidth  - Pixels
%                                Default 1
%                    barLength - Pixels
%                                Default []
%                    gap       - Gap between upper and lower pattern
%                                (filled with bgColor)
%                                Default 0
%                    barColor  - Bar color, e.g., [.6 .4 .2]
%                                Default 1
%                    pattern   - Spatial pattern
%                                Default []
%
%
% Outputs:
%    params   - The vernier stimulus parameters structure
%
% Optional key/value pairs:
%    None.
%
% See Also:
%	 sceneCreate('vernier', ...) , oisCreate('vernier', ...), imageVernier
%

% History:
%    xx/xx/16  BW   ISETBIO Team, 2016
%    02/01/18  jnm  Formatting

% Examples:
%{
    p = vernierP; 
    p.name = 'example';
    s = sceneCreate('vernier', 'display', p);
    ieAddObject(s);
    sceneWindow;

    p.bgColor = [1 0 0];
    p.barColor = [0 1 0]; 
    s = sceneCreate('vernier', 'display', p);
    ieAddObject(s);
    sceneWindow;

    p.barLength = 8;  
    s = sceneCreate('vernier', 'display', p);
    ieAddObject(s);
    sceneWindow;

    x = (-32:32) / 64;
    f = 2;
    pattern = 0.5 * cos(2 * pi * f * x) + 0.5;
    p = vernierP('pattern', pattern, 'offset', 6, 'barLength', 12);
    s = sceneCreate('vernier', 'display', p);
    ieAddObject(s);
    sceneWindow;
%}
%% Parse arguments
p = inputParser;

p.addParameter('name', 'unknown', @ischar);

p.addParameter('display', ...
    displayCreate('LCD-Apple', 'wave', 400:10:700), @isstruct);
p.addParameter('sceneSz', [50 50], @isnumeric);

p.addParameter('offset', 1, @isscalar);
p.addParameter('bgColor', 0.5, @isscalar);

p.addParameter('barWidth', 1, @isnumeric);
p.addParameter('barLength', [], @isscalar);
p.addParameter('barColor', 1, @isscalar);

p.addParameter('gap', 0, @isscalar);
p.addParameter('pattern', [], @ismatrix);

p.parse(varargin{:});

%% Assign parameters
params = p.Results;

% % Identifier
% params.name  = p.Results.name;           % Char
% 
% % Display scene
% params.display  = p.Results.display;     % Display structure
% params.sceneSz  = p.Results.sceneSz;     % Pixels
% 
% % General
% params.offset   = p.Results.offset;      % Pixels
% params.bgColor  = p.Results.bgColor;     % 0 to 1?
% 
% % Bar properties.  Define better
% params.barWidth  = p.Results.barWidth;    % Pixels
% params.barLength = p.Results.barLength;   % Pixels?
% params.barColor  = p.Results.barColor;    % (rgb?)
% 
% params.pattern   = p.Results.pattern;    % Spatial
% params.gap       = p.Results.gap;

end
