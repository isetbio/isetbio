function t_linearFilters
% Compute the photocurrent at different mean field levels
%
% Description:
%   Computes L-, M- and S-cone outer segment photocurrent responses to
%   luminance step stimuli of fixed height presented on different
%   backgrounds. Visualizes isomerization responses, outer segement impulse
%   responses and outersegment photocurrent responses.
%
%   This tutorial mainly shows how the background luminance (adapting
%   stimulus) affects the cone outer-segment linear impulse response
%   (luminance adaptation) and, thus, the ensuing cone photocurrent
%   responses, and how this adaptation depends on cone type.
%
% NOTES
%   * Copying and pasting the cells doesn't work because it has internal
%   functions (e.g., oiGenerate).  We should re-write so the tutorial can
%   be run from the editor.
%
%   * Should we separate the plots by cone type? Should we use
%   vcNewGraphWin? 
%
% NPC, ISETBIO TEAM, 2016
%
% See also: t_osTimeStep

% 01/07/18  npc  Cleaned up, comments

% Define the time axis for the simulation
stimulusSamplingInterval = 10/1000;
oiTimeAxis = -0.5:stimulusSamplingInterval:1;  

% Define the bsckground luminace levels examined and the pedestal step (cd/m2)
backgroundLuminances = [25 50 100 200 400]; 
lumStep = 20;

% Generate optics
noOptics = false;
theOI = oiGenerate(noOptics);

% Generate a cone mosaic with 1 L-, 1 M-, and 1 S-cone only
mosaicSize =  nan;                    
integrationTime = 1/1000;            
photonNoise = 'none';                                        
osNoise = 'none'; 
osTimeStep = 0.2/1000;
theConeMosaic = coneMosaicGenerate(mosaicSize, photonNoise, osNoise, integrationTime, osTimeStep);

modulationFunction = cell(numel(backgroundLuminances),1);
photocurrents      = cell(numel(backgroundLuminances),1);
isomerizations     = cell(numel(backgroundLuminances),1);
osLinearFilters    = cell(numel(backgroundLuminances),1);
for iLum = 1:numel(backgroundLuminances)
    fprintf('Computing os linear filters for background luminance %2.1f cd/m2  [%d/%d]\n',  backgroundLuminances(iLum), iLum, numel(backgroundLuminances));
    
    % Compute scene
    FOV = 1;
    theScene = uniformFieldSceneCreate(FOV, backgroundLuminances(iLum));
    
    % Compute the optical images
    oiBackground = oiCompute(theOI, theScene);
    oiModulated  = oiBackground;

    % Generate the sequence of optical images
    % Step stimulus modulation function
    modulationFunction{iLum} = zeros(1,numel(oiTimeAxis));
    modulationFunction{iLum}(oiTimeAxis>0.2 & oiTimeAxis<0.6) = lumStep/backgroundLuminances(iLum);
    theOIsequence = oiSequence(oiBackground, oiModulated, oiTimeAxis, modulationFunction{iLum});

    % Generate eye movement path with zero movement
    eyeMovementsNum = theOIsequence.maxEyeMovementsNumGivenIntegrationTime(theConeMosaic.integrationTime);
    instancesNum = 1;
    emPaths = zeros(instancesNum, eyeMovementsNum,2);

    % Compute responses
    [isomerizations{iLum}, photocurrents{iLum}, osLinearFilters{iLum}] = ...
            theConeMosaic.computeForOISequence(theOIsequence, ...
            'emPaths', emPaths, ...
            'currentFlag', true);
    timeAxis = theConeMosaic.timeAxis + theOIsequence.timeAxis(1);
end

visualizedTrange = [0 1];
% Visualize the luminance time course
render(backgroundLuminances, modulationFunction, theConeMosaic.integrationTime, 'stimulus luminance (cd/m2)', [0 500], oiTimeAxis, visualizedTrange,  1);

% Visualize the isomerizations time courses
render(backgroundLuminances, isomerizations, theConeMosaic.integrationTime, sprintf('isomerizations/%d ms',theConeMosaic.integrationTime*1000), [0 40], timeAxis, visualizedTrange,  2);

% Visualize the outer-segment linear filters
render(backgroundLuminances, osLinearFilters, theConeMosaic.integrationTime, '', [-0.1 1.0]*0.18, [], [], 3);

% Visualize the outer-segement photocurrent responses
render(backgroundLuminances, photocurrents, theConeMosaic.integrationTime, 'photocurrent (pAmps)', [-90 -40], timeAxis, visualizedTrange, 4);
end

function render(backgroundLuminances, responses, integrationTime, yLabelTitle, yAxisRange, timeAxis, tRange, figNo)    
subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 1, ...
       'colsNum', 3, ...
       'heightMargin',   0.04, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.05, ...
       'rightMargin',    0.00, ...
       'bottomMargin',   0.08, ...
       'topMargin',      0.04);

coneNames = {'L', 'M', 'S'};
colorIR = (jet(numel(backgroundLuminances)));

for iLum = 1:numel(backgroundLuminances)  
    if (iLum == 1)
        legends = {};
        hFig = figure(figNo); clf;
        set(hFig, 'Position', [10 10 1380 650+50*iLum], 'Color', [0.1 0.1 0.1]);
    end

    legends{numel(legends)+1} = sprintf('lum: %2.1f cd/m2', backgroundLuminances(iLum)); %#ok<AGROW>
    color = squeeze(colorIR(iLum,:));
    if (isempty(timeAxis))
        theTimeAxis = (1:length(responses{iLum}))*integrationTime;
    else
        theTimeAxis = timeAxis;
    end
    for coneIndex = 1:3
        subplot('Position', subplotPosVectors(1,coneIndex).v);

        if (~strcmp(yLabelTitle, 'stimulus luminance (cd/m2)'))
            coneResponses = squeeze(responses{iLum});
            if (size(coneResponses,1) == 3)
                response = coneResponses(coneIndex,:);
            else
                response = coneResponses(:,coneIndex);
            end
        else
           response = backgroundLuminances(iLum) * (1+squeeze(responses{iLum}));
        end

        plot(theTimeAxis, response, 'k-', 'Color', color, 'LineWidth', 1.5);
        hold on;

        if isempty(tRange)
            set(gca, 'XLim', [theTimeAxis(1) theTimeAxis(end)])
        else
            set(gca, 'XLim', tRange)
        end

        set(gca, 'FontSize', 14, 'Color', [0.3 0.3 0.3], 'XColor', [0.9 0.9 0.9], 'YColor', [0.9 0.9 0.9]);
        if (~isempty(yAxisRange))
            set(gca, 'YLim', yAxisRange);
        end

        if (coneIndex == 1)
             ylabel(yLabelTitle, 'FontWeight', 'bold');
        end
        xlabel('time (seconds)', 'FontWeight', 'bold');

        grid on; box on;
        hL = legend(legends);
        set(hL, 'FontSize', 14, 'Color', [0.6 0.6 0.6]);
        title(sprintf('%s-cone', coneNames{coneIndex}), 'FontSize', 14, 'FontWeight', 'bold', 'Color', [1 1 1]);
    end % coneIndex
    drawnow;
end %iLum
end

% Internal utility
function theConeMosaic = coneMosaicGenerate(mosaicSize, photonNoise, osNoise, integrationTime, osTimeStep)
% Default human mosaic
theConeMosaic = coneMosaic;

% Adjust size
if isnan(mosaicSize)
    % Generate a human cone mosaic with 1L, 1M and 1S cone
    theConeMosaic.rows = 1;
    theConeMosaic.cols = 3;
    theConeMosaic.pattern = [2 3 4];
else
    theConeMosaic.setSizeToFOV(mosaicSize);
end

% Set the noise
theConeMosaic.noiseFlag = photonNoise;

% Set the integrationTime
theConeMosaic.integrationTime = integrationTime;

% Generate the outer-segment object to be used by the coneMosaic
theOuterSegment = osLinear();
theOuterSegment.noiseFlag = osNoise;

% Set a custom timeStep, for @osLinear we do not need the default 0.1 msec
theOuterSegment.timeStep = osTimeStep;

% Couple the outersegment object to the cone mosaic object
theConeMosaic.os = theOuterSegment;
end


function theOI = oiGenerate(noOptics)
% Generate optics
if (noOptics)
    theOI = oiCreate('diffraction limited');
    optics = oiGet(theOI,'optics');           
    optics = opticsSet(optics,'fnumber',0);
    optics = opticsSet(optics, 'off axis method', 'skip');
    theOI = oiSet(theOI,'optics', optics);
else
    theOI = oiCreate('human');
end
end

% Internal utility.
function uniformScene = uniformFieldSceneCreate(FOV, meanLuminance)
uniformScene = sceneCreate('uniform equal photon', 128);
% square scene with desired FOV
uniformScene = sceneSet(uniformScene, 'wAngular', FOV);
% 1 meter away
uniformScene = sceneSet(uniformScene, 'distance', 1.0);
% adjust radiance according to desired  mean luminance
uniformScene = sceneAdjustLuminance(uniformScene, meanLuminance);
end