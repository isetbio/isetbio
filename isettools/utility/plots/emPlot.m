function emPlot(fixEM,plotType,varargin)
% Fixational eye movement plotting gateway routine
%
% Description
%   We will collect up various methods for plotting the eye movements and
%   visualizing their statistics here.
%
% Inputs:
%   fixEM:
%   plotType: 'emMosaic'
%   
% Returns:
%
% Optional key/value pairs
%   'cone mosaic'
%   'visualized trial'
%
% BW, May 13, 2018

% History:
%   06/21/18  dhb  Fix example, small typo.

% Examples:
%{
% EM superimposed on cone mosaic

sparams.fov = .2;
stimWeights = zeros(1, 50); stimWeights(2:4) = 1;
impulse = oisCreate('impulse', 'add', stimWeights, ...
    'sceneParameters', sparams);

fovDegs = sparams.fov* 0.8;
cm = coneMosaicHex(4,'fovDegs', fovDegs);
cm.integrationTime = 1e-3;   % Secs;
eyeMovementsPerTrial = impulse.maxEyeMovementsNumGivenIntegrationTime(cm.integrationTime);

fixEM = fixationalEM();
fixEM.computeForConeMosaic(cm, eyeMovementsPerTrial, ...
    'nTrials', 1, ...
    'rSeed', 1);
emPlot(fixEM,'emMosaic','conemosaic',cm);
%}

%%
p = inputParser;
varargin = ieParamFormat(varargin);
p.addRequired('fixEM',@(x)(isa(x,'fixationalEM')));
p.addRequired('plotType',@ischar);

p.addParameter('conemosaic',[],@(x)(isa(x,'coneMosaicHex')));
p.addParameter('visualizedtrial',[],@isscalar);

p.parse(fixEM,plotType,varargin{:});

switch(p.Results.plotType)
    case 'emMosaic'
        cm = p.Results.conemosaic;
        visualizedTrial = p.Results.visualizedtrial;
        
        hFig = vcNewGraphWin;
        set(hFig, 'Position', [0 0 0.3 0.5]);

        % The @fixationalEM class can provide the generated emPath in units
        % of 'arc min', 'microns', or 'cone mosaic pattern steps'
        % The @coneMosaicHex.visualizeGrid() method expects emPosMicrons
        cm.visualizeGrid(...
            'axesHandle', gca, ...
            'overlayEMpathmicrons', squeeze(fixEM.emPosMicrons(visualizedTrial,:,:)), ...
            'overlayNullSensors', false, ...
            'apertureShape', 'disks', ...
            'visualizedConeAperture', 'lightCollectingArea', ...
            'labelConeTypes', true,...
            'generateNewFigure', false);
        
        xyRangeMeters = max(...
            [max(max(abs(squeeze(fixEM.emPosMicrons(visualizedTrial,:,:)))))*1e-6 ...
            max(0.5*cm.fov*cm.micronsPerDegree*1e-6)])*1.01;
        
        ticks = (-200:20:200)*1e-6;
        tickLabels = sprintf('%2.0f\n', ticks*1e6);
        
        set(gca, 'XLim', xyRangeMeters*[-1 1], 'YLim', xyRangeMeters*[-1 1], ...
            'XTick', ticks, 'YTick', ticks, ...
            'XTickLabels', tickLabels, 'YTickLabels', tickLabels, ...
            'FontSize', 16);
        
        xlabel('space (microns)'); xlabel('space (microns)'); box on;
        
    otherwise
        error('Unknown plotType %s\n',plotType);
end

end