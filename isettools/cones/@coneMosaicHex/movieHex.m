function uData = movieHex(conemosaicH, varargin)
% Generates a movie of the coneMosaicHex absorptions or current over time.
%
% Syntax:
%   uData = movieHex(conemosaicH, [varargin])
%
% Description:
%    Generates a movie of the coneMosaicHex absorptions or current over
%    time. It is used in coneMosaicHex.window().
%
%    This allows coneMosaicWindow to generate the nice movies from the
%    coneMosaicHex class. This is called from @coneMosaicHex/plot.m which
%    is utilized by coneMosaicWindow.m.
%
% Inputs:
%    conemosaicH - The cone mosaic hex to make a movie from.
%    varargin    - (Optional) Additional information for building the movie
%
% Outputs:
%    uData       - The user data containing the movie.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    t_coneMosaicHex3.m, line 70.
%

% History:
%    10/xx/16  JRG/BW/NPC  (c) isetbio team
%    02/16/18  JNM         Formatting

% Examples:
%{
    % TODO: Get someone to fix example. Is currently returning the error
    % 'Index exceeds matrix dimensions'
    cmosaicH = coneMosaicHex(2);
    movieHex(cmosaicH, 'type', 'current')
%}
%% Parse input data
p = inputParser;
addRequired(p, 'conemosaicH', @(x) isa(x, 'coneMosaicHex'));
addParameter(p, 'type', 'absorptions', @ischar);
p.parse(conemosaicH, varargin{:});
conemosaicH = p.Results.conemosaicH;
plotType = p.Results.type;

%% Select type of data to plot
switch plotType
    case 'absorptions'
        dataHex = conemosaicH.absorptions;
    otherwise
        dataHex = conemosaicH.current;
end

%% Get the nice coneMosaicHex image for average mosaic response over time
% Render activation images for the hex mosaic
tic
disp('Calculating activation density map for hex data.');
[activationsHexImage, ~] = ...
    conemosaicH.computeActivationDensityMap(dataHex);
toc

activationsHexMovie = zeros([size(activationsHexImage), size(dataHex, 3)]);
for frameIndex = 1:size(dataHex, 3)
    [activationsHexImage, ~] = ...
        conemosaicH.computeActivationDensityMap(dataHex(:, :, frameIndex));
    activationsHexMovie(:, :, frameIndex) = activationsHexImage;
end

%% Show movie
% set(gca, 'CLim', isomerizationsRange, 'XTick', [], 'YTick', []);
axis 'image';
axis 'xy';

% title('hex mosaic isomerizations (all cones)', 'FontSize', 16);
% % % % % % %
uData = ieMovie(activationsHexMovie, varargin{:});

%%