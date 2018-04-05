function plotHexMeanImage(conemosaicH, varargin)
% Plot mean absorptions or current from cone mosaic hex into a window image
%
% Syntax:
%    plotHexMeanImage(coneMosaicH, [varargin])
%
% Description:
%    Plots the mean absorptions or current from a coneMosaicHex as a nice
%    image. Used in coneMosaicHex.window().
% 
%    This allows coneMosaicWindow to generate the nice images from the
%    coneMosaicHex class. This is called from @coneMosaicHex/plot.m which
%    is utilized by coneMosaicWindow.m.
%
%    There are examples contained in the code. To access, type 'edit
%    plotHexMeanImage.m' into the Command Window.
%
% Inputs:
%    coneMosaicH - The cone mosaic hex object
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    The pairs are containined in varargin, and can include keys such as:
%
% See Also:
%    t_coneMosaicHex3.m, line 70. (under development)
% 

% History:
%    10/xx/16  JRG/BW/NPC  (c) isetbio team
%    02/16/18  jnm         Formatting

% Examples:
%{
    % Example doesn't work. Error message states 'index exceeds matrix
    % dimensions' inside of computeActivationDensityMap
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
    plotHexMeanImage(cm, 'activationTimeSeries', isomerizations);
%}

%% Parse input
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
else
    switch plotType
        case 'absorptions'
            dataHex = conemosaicH.absorptions;
        otherwise
            dataHex = conemosaicH.current;
    end
end
conesNum = size(dataHex,1);
timeBins = size(dataHex,2);

for frameIndex = 1:timeBins
    dataHexFrame = squeeze(dataHex(:, frameIndex));
    [activationsHexImage, activationImageLMScone, supX, supY] = ...
        conemosaicH.computeActivationDensityMap(dataHexFrame);
    if (frameIndex == 1)
        activationsHexMovie = zeros(size(activationsHexImage,1), size(activationsHexImage,2), timeBins);
    end
    activationsHexMovie(:, :, frameIndex) = activationsHexImage;
end

%% Get the average mosaic response over time
activationsHexImage = squeeze(mean(activationsHexMovie,3));

%% Display results
figure()
imagesc(supX, supY, activationsHexImage);
axis 'image';
axis 'xy';
colormap gray;

hold on;
sampledHexMosaicXaxis = conemosaicH.patternSupport(1, :, 1) + ...
    conemosaicH.center(1);
sampledHexMosaicYaxis = conemosaicH.patternSupport(:, 1, 2) + ...
    conemosaicH.center(2);
dx = conemosaicH.pigment.pdWidth;

axis 'equal';
axis 'xy'
xTicks = [sampledHexMosaicXaxis(1) conemosaicH.center(1) ...
    sampledHexMosaicXaxis(end)];
yTicks = [sampledHexMosaicYaxis(1) conemosaicH.center(2) ...
    sampledHexMosaicYaxis(end)];
xTickLabels = sprintf('%2.0f um\n', xTicks * 1e6);
yTickLabels = sprintf('%2.0f um\n', yTicks * 1e6);
set(gca, 'XTick', xTicks, 'YTick', yTicks, 'XTickLabel', xTickLabels, ...
    'YTickLabel', yTickLabels);
set(gca, 'FontSize', 16, 'XColor', [0.1 0.2 0.9], ...
    'YColor', [0.1 0.2 0.9], 'LineWidth', 1.0);
box on;
grid off;
dataRange = prctile(activationsHexImage(:), [5 95]);

set(gca, 'CLim', dataRange);
set(gca, 'XLim', [sampledHexMosaicXaxis(1) - dx ...
    sampledHexMosaicXaxis(end) + dx]);
set(gca, 'YLim', [sampledHexMosaicYaxis(1) - dx ...
    sampledHexMosaicYaxis(end) + dx]);

drawnow
