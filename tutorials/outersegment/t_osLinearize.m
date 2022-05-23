% Show how small perturbation linear os impulse responses are obtained
%
% Description:
%    Illustrates foveal and peripheral cone outer segment temporal impulse
%    responses using the osBioPhys object at different mean intensities.
%
%    This has a some overlap with t_soFoveaPeriphery. The difference is
%    that this calls through the os linearFilters methods, whereas
%    t_osFoveaPeriphery exposes the underlying calls to the osBioPhys
%    object. In t_osFoveaPeriphery the raw response to an incrmenet in the
%    dark is shown. Here the impuluse response starts at zero and shows the
%    differential response to a stimulus on a steady background.
%
%    Note the different current scales for foveal and peripheral
%    responeses. The foveal response is larger and slower.
%
% See Also:
%    t_osFoveaPeriphery
%

% History:
%    11/XX/16  JRG/BW  (c) Isetbio team
%    08/06/17  dhb     Cleaning pass.

%% Initialize
ieInit;

%% Create a small cone mosaic of L cones
cMosaic = coneMosaic('os', osLinear, 'pattern', [2 2 2]);

%% Setup background levels
% Mean of isomerization stimulus in R*/sec (scaled by time step to be
% placed in bins of photons). The mean rate affects the magnitude of the
% impulse response due to adaptation.
%
% For this tutorial, background rates are the same for all classes of cone.
% Results for the L cones are plotted.
meanIsoArray = [100 500 1000 2000 5000 10000] * cMosaic.integrationTime;

%% Allocate space for the impulse responses
% A litte ugly, we get an impulse response and find its size.
nTimeSamples = size(cMosaic.os.linearFilters(cMosaic), 1);
fovea = zeros(nTimeSamples, length(meanIsoArray));
periphery = zeros(nTimeSamples, length(meanIsoArray));

%%  Loop on different background rates and and compute
os = cMosaic.os;
for ii = 1:length(meanIsoArray)
    % Set mean background for start of current output
    cMosaic.absorptions = repmat(meanIsoArray(ii), 1, 3);

    % Compute outer segment currents for the fovea. Do this by pulling out
    % the first column of return, corresponding to the L cones.
    tmp = os.linearFilters(cMosaic, 'eccentricity', 0);
    fovea(:, ii) = tmp(:, 1);

    % Do it for the periphery
    tmp = os.linearFilters(cMosaic, 'eccentricity', 15);
    periphery(:, ii) = tmp(:, 1);

end

%% Plot
vcNewGraphWin([], 'wide');
timeAxis = os.timeAxis;

% Fovea
subplot(1, 2, 1);
hold on;
plot(timeAxis, fovea, 'LineWidth', 3);
grid on;
xlabel('Time (sec)', 'FontSize', 14);
ylabel('pA', 'FontSize', 14);
set(gca, 'fontsize', 14);
title(sprintf('Foveal IR (R*/sec)'));

% Periphery
subplot(1, 2, 2);
hold on;
plot(timeAxis, periphery, 'LineWidth', 3);
grid on;
xlabel('Time (sec)', 'FontSize', 14);
ylabel('pA', 'FontSize', 14);
set(gca, 'fontsize', 14);
title(sprintf('Peripheral IR (R*/sec)'));

% Add a legend
legend('100', '500', '1000', '2000', '5000', '10000')
