function d = generateCustomDisplay(varargin)
% Generate a custom display
%
% Syntax:
%   d = generateCustomDisplay
%
% Description:
%    Generate a custom display, with optional parameters for setting the
%    - dots per inch
%    - viewing distance
%    - wavelength support
%    - gun SPD
%    - ambient SPD
%    - gammaTable
%
% Inputs:
%    None.
%
% Outputs:
%    d   - the generate display object
%
% Optional key/value pairs:
%    'name'                                                  - String. The name of the display
%    'dotsPerInch'                                           - Scalar. The DPI
%    'viewingDistanceMeters'                                 - Scalar. The viewing distance in meters
%    'wavelengthSupportNanoMeters'                           - [nWaves x 1] matrix of wavelengths
%    'spectralPowerDistributionWattsPerSteradianM2NanoMeter' - [nWaves x 3] matrix of the RGB guns SPDs
%    'ambientSPDWattsPerSteradianM2NanoMeter'                - [nWaves x 1] matrix of the ambient SPD
%    'gammaTable'                                            - [mValues x 3] matrix of LUTs
%    'plotCharacteristics'                                   - Flag indicating whether to plot the display characteristics

% History:
%    12/01/21  npc  Wrote it.

    % Generate the default display
    defaultDisplay = displayCreate;
    
    % Parse input
    p = inputParser;
    p.addParameter('name', 'custom display', @ischar);
    p.addParameter('dotsPerInch', displayGet(defaultDisplay, 'dpi'), @isscalar);
    p.addParameter('viewingDistanceMeters', 0.57, @isscalar);
    p.addParameter('wavelengthSupportNanoMeters', displayGet(defaultDisplay, 'wave'), @(x)(isnumeric(x)&&(size(x,2)==1)));
    p.addParameter('spectralPowerDistributionWattsPerSteradianM2NanoMeter', displayGet(defaultDisplay, 'spd'), @(x)(isnumeric(x)&&(size(x,2)==3)));
    p.addParameter('ambientSPDWattsPerSteradianM2NanoMeter', displayGet(defaultDisplay, 'ambient spd'), @(x)(isnumeric(x)&&(size(x,2)==1)));
    p.addParameter('gammaTable', repmat((linspace(0,1,1024)'), [1 3]), @(x)(isnumeric(x)&&(size(x,2)==3)));
    p.addParameter('plotCharacteristics', false, @islogical);
    p.parse(varargin{:});
    
    % Make sure spectral dimensions match
    assert(size(p.Results.wavelengthSupportNanoMeters,1) == size(p.Results.spectralPowerDistributionWattsPerSteradianM2NanoMeter,1), ...
        'Rows of ''wavelengthSupportNanoMeters'' must match the rows of ''spectralPowerDistributionWattsPerSteradianM2NanoMeter''.');
    
    assert(size(p.Results.wavelengthSupportNanoMeters,1) == size(p.Results.ambientSPDWattsPerSteradianM2NanoMeter,1), ...
        'Rows of ''wavelengthSupportNanoMeters'' must match the rows of ''ambientSPDWattsPerSteradianM2NanoMeter''.');
    
    % Generate display struct
    d = defaultDisplay;
    d = displaySet(d, 'name', p.Results.name);
    d = displaySet(d, 'dpi', p.Results.dotsPerInch);
    d = displaySet(d, 'viewing distance', p.Results.viewingDistanceMeters);
    d = displaySet(d, 'spd', p.Results.spectralPowerDistributionWattsPerSteradianM2NanoMeter);
    d = displaySet(d, 'gamma', p.Results.gammaTable);
    d = displaySet(d, 'ambient spd', p.Results.ambientSPDWattsPerSteradianM2NanoMeter);
    
    if (p.Results.plotCharacteristics)
        w = displayGet(d, 'wave');
        spd = displayGet(d, 'spd');
        ambientSPD = displayGet(d, 'ambient spd');
        hFig = figure();
        set(hFig, 'Position', [10 10 1500 400], 'Color', [1 1 1]);
        subplot(1,4,1);
        stairs(w, spd(:,1)*1e3, 'r-', 'LineWidth', 1.5); hold on;
        stairs(w, spd(:,2)*1e3, 'g-','LineWidth', 1.5);
        stairs(w, spd(:,3)*1e3, 'b-','LineWidth', 1.5);
        stairs(w, ambientSPD*1e3, 'k--', 'LineWidth', 1.5);
        legend({'R', 'G', 'B', 'ambient'})
        xlabel('wavelength (nm)');
        ylabel('power (milliWatts/Sr/m2/nm)');
        title(sprintf('peak luminance: %2.1f cd/m2\ndark luminance: %2.3f cd/m2', displayGet(d, 'peak luminance'), displayGet(d, 'dark luminance')));
        axis 'square';
        set(gca, 'XTick', 300:50:900, 'YTick', 0:1:10,'FontSize', 14);
        grid on;
        box on;
        
        colors = [1 0 0; 0 1 0; 0 0 1];
        for k = 1:3
            subplot(1,4,1+k);
            theGamma = displayGet(d, 'gamma');
            plot((1:size(theGamma,1))/size(theGamma,1), theGamma(:,k), 'r-', 'Color', colors(k,:), 'LineWidth', 1.5);
            set(gca, 'XLim', [0 1], 'XTick', 0:0.1:1, 'YLim', [0 1], 'YTick', 0:0.2:1);
            axis 'square';
            grid on
        end
    end
    
end