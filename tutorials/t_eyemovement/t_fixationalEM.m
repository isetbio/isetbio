%% t_fixationalEM
% Show how to generate fixational eye movements and how to do this 
% for a cone mosaic.
%
% Description:
%    * Calculate em straight or using a cone mosaic
%
% Notes:
%

% History:
%    05/16/18  BW   ISETBIO Team, May 16 2018
%    07/16/18  jnm  Formatting

%% Initialize
em = fixationalEM;              % Instantiate a fixationalEM object
em.microSaccadeType = 'none';   % No microsaccades, just drift

%%
duration = 0.5;   % 0.5 second duration
tStep = 0.001;    % time step: 1 msec
nTrials = 50;
computeVelocity = false;
em.compute(duration, tStep, nTrials, computeVelocity);

%%
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

%%
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
