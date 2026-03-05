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
%    If you are starting with a PTB calibration file to get the data, note
%    that the PTB power units convention is
%    WattsPerSteradianM2WavelengthBand, so that you should divide by the
%    wavelength spacing if you're starting with spectra in a PTB
%    calibration file, or as measured by PTB routine MeasSpd. Here we are
%    starting with an ISETBio default display so that we are inside the
%    ISETBio world where power is specified per nanometer rather than per
%    wavelength band.
%
% Inputs:
%    None.
%
% Outputs:
%    d   - The generated display object
%
% Optional key/value pairs:
%    'name'                                                  - String. The name of the display
%    'dotsPerInch'                                           - Scalar. The DPI
%    'viewingDistanceMeters'                                 - Scalar. The viewing distance in
%                                                              meters. Default 0.57.
%    'wavelengthSupportNanoMeters'                           - [nWaves x 1] matrix of wavelengths
%    'spectralPowerDistributionWattsPerSteradianM2NanoMeter' - [nWaves x 3] matrix of the RGB guns SPDs
%    'ambientSPDWattsPerSteradianM2NanoMeter'                - [nWaves x 1] matrix of the ambient SPD
%    'gammaTable'                                            - [mValues x 3] matrix of LUTs
%    'plotCharacteristics'                                   - Flag indicating whether to plot the display
%                                                              characteristics. Default false.
%
% See also: ptb.GenerateIsetbioDisplayObjectFromPTBCalStruct,
%           ptb.GeneratePTCalStructFromIsetbioDisplayObject

% History:
%    12/01/21  npc  Wrote it.
%    01/16/22  dhb  More comments.
%    03/10/22  npc  Now setting custom wavelength support

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
    p.addParameter('visualizationAxes', {}, @(x)(isempty(x)||(iscell(x))));
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
    d = displaySet(d, 'wave', p.Results.wavelengthSupportNanoMeters);
    d = displaySet(d, 'spd', p.Results.spectralPowerDistributionWattsPerSteradianM2NanoMeter);
    d = displaySet(d, 'gamma', p.Results.gammaTable);
    d = displaySet(d, 'ambient spd', p.Results.ambientSPDWattsPerSteradianM2NanoMeter);
    
    
   % lineColors = brewermap(10, 'spectral');

    lineColors = [...
        1 0.3 0.4;
        0.2 0.8 0.4;
        0.4 0.5 1.0];

    if (p.Results.plotCharacteristics)

        w = displayGet(d, 'wave');
        spd = displayGet(d, 'spd');
        ambientSPD = displayGet(d, 'ambient spd');

        if (isempty(p.Results.visualizationAxes))
            hFig = figure();
            set(hFig, 'Position', [10 10 1500 400], 'Color', [1 1 1]);
            axSPD = subplot(1,4,1);
        else
            axSPD = p.Results.visualizationAxes{1};
        end

        cla(axSPD);

        faceColor = lineColors(1,:);
        makeShadedPlot(w, spd(:,1)*1e3, faceColor, faceColor*0.5, axSPD);

        
        hold(axSPD, 'on');
        faceColor = lineColors(2,:);
        makeShadedPlot(w, spd(:,2)*1e3, faceColor, faceColor*0.5, axSPD);
        faceColor = lineColors(3,:);
        makeShadedPlot(w, spd(:,3)*1e3, faceColor, faceColor*0.5, axSPD);
        stairs(axSPD,w, ambientSPD*1e3, 'k--', 'LineWidth', 1.5);
        maxSPD  = max(spd(:)*1e3);
        if (maxSPD < 0.1)
            yTicks = 0:0.01:maxSPD;
        elseif (maxSPD < 0.5)
            yTicks = 0:0.05:maxSPD;
        elseif (maxSPD < 1)
            yTicks = 0:0.1:maxSPD;
        elseif (maxSPD < 5)
            yTicks = 0:0.5:maxSPD;
        elseif (maxSPD < 10)
            yTicks = 0:1:maxSPD;
        else
            yTicks = 0:10:maxSPD;
        end
        legend(axSPD,{'R', 'G', 'B', 'ambient'})
        xlabel(axSPD,'wavelength (nm)');
        ylabel(axSPD,'power (milliWatts/Sr/m2/nm)');
        title(axSPD,...
            sprintf('peak luminance: %2.1f cd/m2\ndark luminance: %2.2f cd/m2', displayGet(d, 'peak luminance'), displayGet(d, 'dark luminance')), ...
            'FontWeight', 'Normal');
        set(axSPD,'XTick', 300:50:900, 'YTick', yTicks,'FontSize', 14);
        set(axSPD, 'XLim', [w(1) w(end)]);
        grid(axSPD,'on');
        box(axSPD, 'on');
        xtickangle(axSPD,0);
        hold(axSPD, 'off');

        colors = [1 0 0; 0 1 0; 0 0 1];
        for k = 1:3
            if (isempty(p.Results.visualizationAxes))
                axRGBgunLUT = subplot(1,4,1+k);
            else
                axRGBgunLUT = p.Results.visualizationAxes{k+1};
            end
            if (~isempty(axRGBgunLUT))
                theGamma = displayGet(d, 'gamma');
                plot(axRGBgunLUT, (1:size(theGamma,1))/size(theGamma,1), theGamma(:,k), 'r-', 'Color', colors(k,:), 'LineWidth', 1.5);
                set(axRGBgunLUT, 'XLim', [0 1], 'XTick', 0:0.2:1, 'YLim', [0 1], 'YTick', 0:0.2:1, 'FontSize', 14);
                xlabel(axRGBgunLUT, 'settings value');
                if (k == 1) || (~isempty(p.Results.visualizationAxes))
                    ylabel(axRGBgunLUT, 'primary value');
                end
                axis(axRGBgunLUT, 'square');
                grid(axRGBgunLUT, 'on');
                xtickangle(axRGBgunLUT, 0);
            end

        end
    end
    
end


function makeShadedPlot(x,y, faceColor, edgeColor, ax)
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    px = [px(1) px px(end)];
    py = [1*eps py 2*eps];
    pz = -10*eps*ones(size(py)); 
    patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor',edgeColor, 'FaceAlpha', 0.5);
end
