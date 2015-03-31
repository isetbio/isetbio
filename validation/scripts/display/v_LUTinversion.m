function varargout = v_LUTinversion(varargin)
%
% Validate display calibration and compare against PTB answers.
%

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

    %% Initialize ISET
    s_initISET;
    
    %% Remove the Brainard lab PTB overrides folder from the path
    [removedFolderFromCurrentPath, originalPath] = removeBrainardLabPTBOverrides();
    
    try
        
        %% Compare PTB vs ISETBIO LUT inversion for three different displays
        displaysToTest = {'OLED-Sony', 'LCD-Apple', 'CRT-Dell'};
        
        %% and four different gamma table lengths
        gammaTableLengthsToTest = [128 256 1024 2048];
        
        for displayIndex = 1:numel(displaysToTest)
            
            %% Create a display object
            d = displayCreate(displaysToTest{displayIndex});

            %% Retrieve key properties of the display
            % Gamma table
            gammaTable = displayGet(d, 'gamma table');
            % Remove 4-th primary, if it exists
            if (size(gammaTable,2) > 3)
                gammaTable = gammaTable(:,1:3);
            end
            originalGammaTableLength = size(gammaTable,1);
            originalSettingsValues   = linspace(0,1,originalGammaTableLength);
            
            % Primary SPD
            wave = displayGet(d, 'wave');
            spd  = displayGet(d, 'spd primaries');
            % Remove 4-th primary, if it exists
            if (size(spd ,2) > 3)
                spd = spd(:,1:3);
            end
            
            % Screen size and DPI
            dotsPerMeter = displayGet(d, 'dots per meter');
            screenSizeInPixels = [1920 1080];
                
            %% Generate PTB-compatible calStruct describing the display
            PTBcal = ptb.GeneratePsychToolboxCalStruct(...
                'name', displayGet(d, 'name'), ...
                'gammaTable', gammaTable, ...
                'wave', wave, ...
                'spd', spd, ...
                'screenSizeInPixels', screenSizeInPixels, ...
                'dotsPerMeter', dotsPerMeter ...
            );
                
            for resolutionIndex = 1:numel(gammaTableLengthsToTest)
                %% Compute the inverse gamma table with desired resolution
                nInputLevels = gammaTableLengthsToTest(resolutionIndex);
                
                % PTB-solution
                PTBcal = CalibrateFitGamma(PTBcal, nInputLevels);
                gammaMethod = 1;
                PTBcal = SetGammaMethod(PTBcal, gammaMethod, nInputLevels);
                PTBinverseGamma       = PTBcal.iGammaTable;
                settingsValues        = PTBcal.gammaInput;

                % ISETBIO solution
                ISETBIOinverseGamma = displayGet(d, 'inverse gamma', nInputLevels);
                % ISETBIOinverseGamma = ieLUTInvert(gammaTable,nInputLevels);
                ISETBIOinverseGamma = ISETBIOinverseGamma / (originalGammaTableLength-1);
                
                % Log the data 
                dataStruct = struct( ...
                    'displayName', displaysToTest{displayIndex}, ...
                    'originalSettingsValues', originalSettingsValues, ...
                    'gammaTable', gammaTable, ...
                    'nInputLevels', nInputLevels, ...
                	'settingsValues', settingsValues, ...
                	'PTBinverseGamma', PTBinverseGamma, ...
                	'ISETBIOinverseGamma', ISETBIOinverseGamma ...
                 );
                
                data{displayIndex}{resolutionIndex} = dataStruct;
            end
        end
        
    catch err
        % restore original path
        if (~isempty(removedFolderFromCurrentPath))
            path(originalPath);
        end
        
        % retrhow the error
        rethrow(err);
    end
    
    % restore original path
    if (~isempty(removedFolderFromCurrentPath))
       path(originalPath);
    end
    
    %% Plot
    if (runTimeParams.generatePlots)
        h = figure(1);
        set(h, 'Position', [10 10 1704 1196]);
        clf;
        
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
    end
end

% Helper method to remove the BrainardLabPTBOverrides folder if it exists on the current path.
% We want to remove this override, so we can use the original PTB functions
% without the CalStructOBJ (which is found only on our computers).
function [removedFolderFromCurrentPath, originalPath] = removeBrainardLabPTBOverrides

    removedFolderFromCurrentPath = '';
    originalPath = path;
    
    % Folder to remove the Overrides/PTB-3 folder from the current path, if the folder exists on the path
    PTBoverridesDirToRemoveFromPath = '/Users/Shared/Matlab/Toolboxes/BrainardLabToolbox/Overrides/PTB-3';
    
    % determine if the PTBoverridesDirToRemoveFromPath is in the current path
    pathCell = regexp(path, pathsep, 'split');
    onPath = any(strcmpi(PTBoverridesDirToRemoveFromPath, pathCell));
    if (onPath)
        rmpath(PTBoverridesDirToRemoveFromPath);
        removedFolderFromCurrentPath = PTBoverridesDirToRemoveFromPath;
    end
end

