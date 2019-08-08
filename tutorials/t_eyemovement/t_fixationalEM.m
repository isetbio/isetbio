%% t_fixationalEM
% This tutorial demonstrates how to generate fixational eye movements. It
% first generates 50 different eye movements sequences, each over the
% specified duration and temporal sampling. These 50 eye movements are
% displayed in a figure. Next, the eye movements are re-generated and
% re-plotted, but this time spatially sampled on a cone mosaic. Lastly, one
% of the eye movement sequences is chosen and used to generate a movie of
% cone excitations from a simple example scene. 
%
% Description:
%    * Calculate em straight or using a cone mosaic
%
% Notes:
%

% History:
%    05/16/18  BW   ISETBIO Team, May 16 2018
%    07/16/18  jnm  Formatting
%    07/29/19  TL   Clarification and improvements.

%% Initialize
em = fixationalEM;              % Instantiate a fixationalEM object
em.microSaccadeType = 'none';   % No microsaccades, just drift

%% Specify parameters and generate an eye movement sequence
duration = 0.25;   % 0.25 second duration
tStep = 0.001;    % time step: 1 msec
nTrials = 50;
computeVelocity = false;
em.compute(duration, tStep, nTrials, computeVelocity);

%% Plot the eye positions
visFOVArcMin = 15;  % visualized area in arc min
vcNewGraphWin;
hold on
for ii = 1:nTrials
    thisPath = squeeze(em.emPosArcMin(ii, :, :));
    plot(thisPath(:, 1), thisPath(:, 2), '.-');
end
grid on;
set(gca, 'XLim', visFOVArcMin*[-1 1], 'YLim', visFOVArcMin*[-1 1]);
axis 'square'
xlabel('Pos (arc min)');
ylabel('Pos (arc min)');
title(sprintf('Duration %.1f', duration));
drawnow;

%% Show how it is done for a cone mosaic
cm = coneMosaic;
cm.integrationTime = tStep;
nEyeMovements = duration / tStep;
em.computeForConeMosaic(cm, nEyeMovements, 'nTrials', nTrials);

%% Plot the eye positions sampled to the cone mosaic grid
vcNewGraphWin;
hold on
conesToArcMin = cm.patternSampleSize*1e6/300*60;
for ii = 1:nTrials
    thisPath = squeeze(em.emPos(ii, :, :));
    plot(thisPath(:, 1)*conesToArcMin, thisPath(:, 2)*conesToArcMin, '.-');
end
grid on;
set(gca, 'XLim', visFOVArcMin*[-1 1], 'YLim', visFOVArcMin*[-1 1]);
axis 'square'
xlabel('Pos (arc min)')
ylabel('Pos (arc min)')
title(sprintf('Duration %.1f', duration));

%% Generate cone isomerization movie 
% Use one of the above generated eye movements.

% Specify the eye movement positions in the cone mosaic object.
cm.emPositions = squeeze(em.emPos(1,:,:));

% Create a scene...
s = sceneCreate('rings rays');
s = sceneSet(s, 'fov', 1);

% ...then the oi
oi = oiCreate;
oi = oiCompute(oi, s);

% Calculate cone isomerizations
cm.compute(oi);

cm.window();

% Show the video either through the cone mosaic drop down menu or with the
% following command:
% cm.plot('movieabsorptions')
