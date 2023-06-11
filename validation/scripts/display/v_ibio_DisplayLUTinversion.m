% Validate display calibration LUT inversion between ISETBio and PTB.
%
% The inversion methods are not identical.  In particular, the Psychtoolbox
% fit the gamma function with a power function, which smooths things out.
% (The PTB functions have many other options for how they fit gamma
% functions, but we are not exploring those here.)  The ISETBio functions
% do not do the smoothing.
%

%% Compare PTB vs ISETBIO LUT inversion for three different displays
% and four different gamma table lengths
displaysToTest = {'OLED-Sony', 'LCD-Apple', 'CRT-Dell'};
gammaTableLengthsToTest = [128 256 1024 2048];
nDisplays = numel(displaysToTest);

% Loop over displays
ieNewGraphWin([],'wide');
for displayIndex = 1:nDisplays

    % Create a display object
    d = displayCreate(displaysToTest{displayIndex});

    % Retrieve key properties of the display
    %
    % Gamma table, remove 4-th primary, if it exists,
    % for testing purposes.
    gammaTable = displayGet(d, 'gamma table');
    if (size(gammaTable,2) > 3)
        gammaTable = gammaTable(:,1:3);
    end
    originalGammaTableLength = size(gammaTable,1);
    originalSettingsValues   = linspace(0,1,originalGammaTableLength);

    % Primary SPDs
    % Remove 4-th primary, if it exists, for testing purposes.
    wave = displayGet(d, 'wave');
    spd  = displayGet(d, 'spd primaries');

    % Screen size and DPI
    dotsPerMeter = displayGet(d, 'dots per meter');
    screenSizeInPixels = [1920 1080];

    % Generate PTB-compatible calStruct describing the display
    PTBcal = ptb.GeneratePsychToolboxCalStruct(...
        'name', displayGet(d, 'name'), ...
        'gammaInput', originalSettingsValues, ...
        'gammaTable', gammaTable, ...
        'wave', wave, ...
        'spd', spd, ...
        'ambientSpd', zeros(length(wave),1), ...
        'screenSizeInPixels', screenSizeInPixels, ...
        'dotsPerMeter', dotsPerMeter ...
        );

    for resolutionIndex = 1:numel(gammaTableLengthsToTest)
        % Compute the inverse gamma table with desired resolution
        nInputLevels = gammaTableLengthsToTest(resolutionIndex);

        % PTB-solution
        PTBcal = CalibrateFitGamma(PTBcal, nInputLevels);
        gammaMethod = 1;
        PTBcal = SetGammaMethod(PTBcal, gammaMethod, nInputLevels);
        PTBinverseGamma       = PTBcal.iGammaTable;
        settingsValues        = PTBcal.gammaInput;

        % ISETBIO solution
        ISETBIOinverseGamma = displayGet(d, 'inverse gamma', nInputLevels);

        % Normalize to max of 1 for comparison with PTB
        ISETBIOinverseGamma = ISETBIOinverseGamma / (originalGammaTableLength-1);
       
        subplot(1,nDisplays,displayIndex)
        plot(PTBinverseGamma(:),ISETBIOinverseGamma(:),'.');
        hold on;
        
        %         % Log the data
        %         dataStruct = struct( ...
        %             'displayName', displaysToTest{displayIndex}, ...
        %             'originalSettingsValues', originalSettingsValues, ...
        %             'gammaTable', gammaTable, ...
        %             'nInputLevels', nInputLevels, ...
        %             'settingsValues', settingsValues, ...
        %             'PTBinverseGamma', PTBinverseGamma, ...
        %             'ISETBIOinverseGamma', ISETBIOinverseGamma ...
        %             );
        %
        %         data{displayIndex}{resolutionIndex} = dataStruct;
        %     end
    end
    xlabel('PTB'); ylabel('ISETBio'); axis equal; grid on; identityLine;
    title(displayGet(d,'name'));
        
end

% The values in the dark areas differ, possibly because of the measurement
% noise that is removed by the polynomial function used by PTB.

%% Plot needs some simplification

%{
h = ieNewGraphWin([],'upper left big');
subplotIndex = 0;
for displayIndex = 1:numel(displaysToTest)
    for resolutionIndex = 1:numel(gammaTableLengthsToTest)

        dataStruct = data{displayIndex}{resolutionIndex};
        subplotIndex = subplotIndex + 1;
        subplot(numel(displaysToTest),numel(gammaTableLengthsToTest),subplotIndex);
        hold on;

        % Plot results for RED channel only
        channelIndex = 1;

        % original gamma
        plot(dataStruct.originalSettingsValues, dataStruct.gammaTable(:,channelIndex), 'g.');

        % PTB-inverted
        plot(dataStruct.settingsValues, dataStruct.PTBinverseGamma(:,channelIndex), 'r.');

        % ISETBIO-inverted
        plot(dataStruct.settingsValues, dataStruct.ISETBIOinverseGamma(:,channelIndex), 'b.');

        % ISETBIO-inverted
        plot(dataStruct.settingsValues, 0.5+dataStruct.PTBinverseGamma(:,channelIndex)-dataStruct.ISETBIOinverseGamma(:,channelIndex), 'k-', 'LineWidth', 2.0);

        hold off;

        h_legend = legend('Orig. LUT', 'PTB inverse', 'ISETBIO inverse', '0.5+(PTB-ISETBIO)', 'Location', 'NorthWest');
        set(h_legend,'FontName', 'Menlo', 'FontSize',12);

        set(gca, 'XLim', [0 1], 'YLim', [0 1], 'XTick', [0:0.1:1.0], 'YTick', [0:0.1:1.0], 'FontSize', 12, 'FontName', 'Helvetica');
        xlabel('gamma in','FontSize', 12, 'FontName', 'Helvetica', 'FontWeight', 'bold');
        ylabel('gamma out', 'FontSize', 12, 'FontName', 'Helvetica', 'FontWeight', 'bold');
        axis 'square'; box on; grid on;
        title(sprintf('%s - (%d levels)', dataStruct.displayName, dataStruct.nInputLevels), 'FontSize', 16, 'FontName', 'Helvetica', 'FontWeight', 'bold');
        hold off
    end
end
%}