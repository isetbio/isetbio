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
%e
% Outputs:
%    uData       - The user data containing the movie.
%
% Optional key/value pairs:
%    'activationTimeSeries'   - specific time series to visualize. if not 
%                               passed, data are taken from the @coneMosaic
%                               object
%    'type'                   - String. Specifies the data source from the
%                               @coneMosaic ('absorptions' or 'photocurrents')
%
% See Also:
%    t_coneMosaicHex3.m, line 70.
%

% History:
%    10/xx/16  JRG/BW/NPC  (c) isetbio team
%    02/16/18  JNM         Formatting
%    04/05/18  NPC         Fixed broken example

% Examples:
%{
    hparams(2) = harmonicP;
    hparams(2).freq = 8;
    hparams(2).GaborFlag = .2;
    hparams(1) = hparams(2);
    hparams(1).contrast = 0;
    sparams.fov = 0.5;
    stimWeights = ieScale(fspecial('gaussian', [1, 50], 15), 0, 1);
    ois = oisCreate('harmonic', 'blend', stimWeights, ...
        'testParameters', hparams, 'sceneParameters', sparams);
    cm = coneMosaicHex(7, 'fovDegs', sparams.fov);
    nTrials = 1;
    nEyeMovements = ois.maxEyeMovementsNumGivenIntegrationTime(cm.integrationTime)
    emPaths = zeros(nTrials, nEyeMovements, 2);
    isomerizations = cm.computeForOISequence(ois, 'emPaths', emPaths);
    movieHex(cm, 'activationTimeSeries', isomerizations);
%}
%% Parse input data
p = inputParser;
addRequired(p, 'conemosaicH', @(x) isa(x, 'coneMosaicHex'));
addParameter(p, 'activationTimeSeries', [], @isnumeric);
addParameter(p, 'type', 'absorptions', @ischar);
p.parse(conemosaicH, varargin{:});
conemosaicH = p.Results.conemosaicH;
plotType = p.Results.type;

%% Select type of data to plot
if (~isempty(p.Results.activationTimeSeries))
    trialToVisualize = 1;
    dataHex = squeeze(p.Results.activationTimeSeries(trialToVisualize,:,:));
    conesNum = size(dataHex,1);
    timeBins = size(dataHex,2);
else
    switch plotType
        case 'absorptions'
            dataHex = conemosaicH.absorptions;
        otherwise
            dataHex = conemosaicH.current;
    end
    activeConeIndices = find(conemosaicH.pattern>1);
    timeBins = size(dataHex,3); 
    dataHex = RGB2XWFormat(dataHex);
    dataHex = dataHex(activeConeIndices,:);
    conesNum = size(dataHex,1);
end


%% Get the nice coneMosaicHex image for average mosaic response over time
% Render activation images for the hex mosaic
tic
disp('Calculating activation density map for hex data.');
timeBin = 1;
[activationsHexImage, ~] = ...
    conemosaicH.computeActivationDensityMap(squeeze(dataHex(:,timeBin)));

activationsHexMovie = zeros([size(activationsHexImage), timeBins]);
for frameIndex = 1:timeBins
    dataHexFrame = squeeze(dataHex(:, frameIndex));
    [activationsHexImage, ~] = ...
        conemosaicH.computeActivationDensityMap(dataHexFrame);
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