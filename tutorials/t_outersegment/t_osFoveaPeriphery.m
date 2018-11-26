% Illustrate difference in photocurrent response between fovea & periphery
%
% Description:
%    Illustrates the difference between the foveal and peripheral cone
%    outer segment responses, using the osBioPhys and the osLinear objects.
%
%    A single cone is created, and its absorption time course is set to
%    have an impulse at the first time step. Biophysical and linear outer
%    segments are created, one segment with foveal dynamics and one with
%    peripheral dynamics.
%
%    This tutorial is a variation of validation routine v_osBioPhysObject.
%    It works by generating a cone mosaic and directly setting the photon
%    absorptions there. Better would be to create a stimulus movie and push
%    the example all the way from the stimulus to the photocurrent.
%
%    There is a method for the outersegment object, linearFilters, that
%    gets the impulse response for a small perturbation around a steady
%    uniform background. Tutorial t_osLinearize shows how to do roughly
%    the same thing as here but using that routine. In that tutorial, what
%    you get is the differential response to a stimulus on the steady
%    background, whereas here the raw current response to a flash in the
%    dark is shown.
%
% See Also:
%    t_osLinearize

% History:
%    10/XX/16  JRG (c) Isetbio team
%    08/05/17  dhb  Clean up
%                   Rewrite to use new key/value interface to osBioPhys
%    08/17/17  npc  Added t_osLinear test.
%    11/19/18  JNM  Formatting

%% Initialize
ieInit;

%% Set up an impluse stimulus
% Stimulus parameters:
%    nSamples:    Number of time samples
%    timeStep:    Time step
%    flashIntens: The flash intensity in R*/cone/sec. This is maintained
%                 for 1 bin only.
nSamples = 2000;
timeStep = 1e-4;
flashIntens = 50000;

%% Create stimulus as a 1 by 1 pixel movie
stimulus = zeros(nSamples, 1);
stimulus(1) = flashIntens * timeStep;
stimulus = reshape(stimulus, [1 1 nSamples]);

%% Loop over both outer segment types
osTypesExamined = {'osLinear', 'osBioPhys'};
figPositions = {'upper left', 'upper right'};
for osIndex = 1:numel(osTypesExamined)
    % Get osType and figPosition
    osType = osTypesExamined{osIndex};
    figPosition = figPositions{osIndex};

    % Generate a peripheral cone mosaic
    %    Create an os object with peripheral parameters and insert it into
    %    a cone mosaic. Passing a 1 by 1 pattern matrix to the mosaic
    %    create call causes the mosaic to have one L cone (because pattern
    %    type 2 -> L cone).
    osPeripheral = eval(sprintf('%s(''eccentricity'', 15)', osType));
    osPeripheral.set('noise flag', 'none');
    cmPeripheral = coneMosaic('os', osPeripheral, 'pattern', 2);
    cmPeripheral.integrationTime = timeStep;
    cmPeripheral.os.timeStep = timeStep;

    % Put in our movie of absorptions created above. As noted in the
    % header comments, it would be cleaner to start with a one pixel
    % movie out in scene space and push this through an oi and then
    % compute absorptions from that.
    cmPeripheral.absorptions  = stimulus;

    % Use the computeCurrent method on the cone mosaic object
    %    Method knows it is producing photocurrent from the isomerizations.
    cmPeripheral.computeCurrent;

    % Plot the peripheral current against time
    vcNewGraphWin([], figPosition);
    hold on
    timeAxis = (1:nSamples) * timeStep;
    plot(timeAxis, squeeze(cmPeripheral.current), 'r', 'LineWidth', 2);

    % Repeat for foveal dynamics and add to plot
    %    Setting 'eccentricity' to 0, for foveal dynamics
    osFoveal = eval(sprintf('%s(''eccentricity'', 0)', osType));
    osFoveal.set('noise flag', 'none');
    cmFoveal = coneMosaic('os', osFoveal, 'pattern', 2); % a single cone
    cmFoveal.integrationTime = timeStep;
    cmFoveal.os.timeStep = timeStep;
    cmFoveal.absorptions = stimulus;
    cmFoveal.computeCurrent();
    current2Scaled = (cmFoveal.current) - cmFoveal.current(1);

    % Add foveal to the plot
    %    [Note that the dynamics are slower in the fovea, but that the
    %    amplitude is bigger.]
    plot(timeAxis, squeeze(cmFoveal.current), 'b', 'LineWidth', 2);
    grid on;
    xlabel('Time (sec)', 'FontSize', 14);
    ylabel('Photocurrent (pa)', 'FontSize', 14);
    title(sprintf("'%s' impulse response in the dark", osType), ...
        'FontSize', 16);
    set(gca, 'fontsize', 14);
    legend('Peripheral', 'Foveal');

end
