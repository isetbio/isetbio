%% t_fixationalEMParameters
%
% * Calculate em straight or using a cone mosaic
% * Effects of changing various fixational EM parameters,
%
% Notes:
%   The random seed doesn't seem to work between compute() and
%   computeForConeMosaic().  Look into it.
%
% BW, ISETBIO Team,  May 16 2018

%% Initialize

em = fixationalEM;
em.microSaccadeType = 'none';
em.randomSeed = 1;

%%
duration = 1;     % Secs
tStep   = 0.001; % Secs
nTrials = 5;
computeVelocity = false;
em.compute(duration,tStep,nTrials,computeVelocity);

%%
vcNewGraphWin;
hold on
for ii=1:nTrials
    thisPath = squeeze(em.emPosArcMin(ii,:,:));
    plot(thisPath(:,1),thisPath(:,2),'-');
end
grid on; axis equal
xlabel('Pos (arc min)');
ylabel('Pos (arc min)');
title(sprintf('Duration %.1f',duration));
drawnow;

%% Show how it is done for a cone mosaic

cm = coneMosaic;
nEyeMovements = duration/tStep;
em.computeForConeMosaic(cm,nEyeMovements, ...
    'nTrials',nTrials, ...
    'rSeed', 1);

%%
vcNewGraphWin;
hold on
for ii=1:nTrials
    thisPath = squeeze(em.emPos(ii,:,:));
    plot(thisPath(:,1),thisPath(:,2),'-');
end
grid on; axis equal
xlabel('Pos (cone)')
ylabel('Pos (cone)')
title(sprintf('Duration %.1f',duration));

%%
