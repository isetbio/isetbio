function uData = movieHex(conemosaicH, varargin)
% Generates a movie of the coneMosaicHex absorptions or current over time. 
% Used in coneMosaicHex.window().
% 
%       movieHex(cmosaicH,'type','current')
% 
% This allows coneMosaicWindow to generate the nice movies from the
% coneMosaicHex class. This is called from @coneMosaicHex/plot.m which is
% utilized by coneMosaicWindow.m.
% 
% See also t_coneMosaicHex3.m, line 70.
% 
% 10/2016 (c) JRG/BW/NPC isetbio team

%% Parse input data
p = inputParser;
addRequired(p,'conemosaicH',@(x) isa(x, 'coneMosaicHex'));
addParameter(p,'type','absorptions',@ischar);
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

%% Get the nice coneMosaicHex image for the average mosaic response over time

% Render activation images for the hex mosaic
tic
disp('Calculating activation density map for hex data.');
[activationsHexImage, ~] = conemosaicH.computeActivationDensityMap(dataHex);
toc

activationsHexMovie = zeros([size(activationsHexImage),size(dataHex,3)]);
for frameIndex = 1:size(dataHex,3)
    [activationsHexImage, ~] = conemosaicH.computeActivationDensityMap(dataHex(:,:,frameIndex));
    activationsHexMovie(:,:,frameIndex) = activationsHexImage;
end

%% Show movie
% set(gca, 'CLim', isomerizationsRange, 'XTick', [], 'YTick', []);
axis 'image'; axis 'xy'

% title('hex mosaic isomerizations (all cones)', 'FontSize', 16);
% % % % % % %
uData = ieMovie(activationsHexMovie,varargin{:});

%%