%%t_eyeMovements  Tutorial on controlling eye movement parameters
%
% Deprecated.  We are changing the eye movement model.



%% Build an eye movement  and cone mosaic structure
em = emCreate;
cm = coneMosaic;

% Generate an eye movement sequence nd plot it
cm.emGenSequence(100,'em',em);
vcNewGraphWin;
plot(cm.emPositions(:, 1), cm.emPositions(:, 2));
grid on; xlabel('Horizontal position (cones)');
ylabel('Vertical position (cones)');
set(gca,'xlim',[-30 30],'ylim',[-30 30]);

%% Increase the tremor amplitude

em = emSet(em,'tremor amplitude', 0.0073*3);
cm.emGenSequence(100,'em',em);
vcNewGraphWin;
plot(cm.emPositions(:, 1), cm.emPositions(:, 2));
grid on; xlabel('Horizontal position (cones)');
ylabel('Vertical position (cones)');
set(gca,'xlim',[-30 30],'ylim',[-30 30]);

%% Reset the tremor amplitude and increase the drift
em = emSet(em,'tremor amplitude', 0.0073);
em = emSet(em,'drift speed', 8.7266e-04 * 5);

cm.emGenSequence(100,'em',em);
vcNewGraphWin;
plot(cm.emPositions(:, 1), cm.emPositions(:, 2));
grid on; xlabel('Horizontal position (cones)');
ylabel('Vertical position (cones)');
set(gca,'xlim',[-30 30],'ylim',[-30 30]);

%%  Large drift and tremor, as if in a bad eye movement system

em = emSet(em,'drift speed', 8.7266e-04 * 5);
em = emSet(em,'tremor amplitude', 0.0073*3);
cm.emGenSequence(100,'em',em);
vcNewGraphWin;
plot(cm.emPositions(:, 1), cm.emPositions(:, 2));
grid on; xlabel('Horizontal position (cones)');
ylabel('Vertical position (cones)');
set(gca,'xlim',[-30 30],'ylim',[-30 30]);

%%
