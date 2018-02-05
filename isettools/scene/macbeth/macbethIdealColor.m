function target = macbethIdealColor(illuminant, colorSpace)
% Calculate MCC values for given an illuminant in a color space
%
% Syntax:
%	target = macbethIdealColor(illuminant, colorSpace)
%
% Description:
%    Calculate the macbeth color chart values for a given illuminant in a
%    specified color space.
%
%    Possible illuminants are listed in the function illuminantRead if you
%    feel they are not adequately covered here.
%
%    Examples of function usage are located in the code. To access, enter
%    'edit macbethIdealColor.m' into the Command Window
%
% Inputs:
%    illuminant - (Optional) The illuminant source. It can be a single
%                 entity, or a structure. Default is 'd65'. Options for
%                 either the singular entity, or the structure's .name
%                 component are:
%                    '555nm'          Monochrome
%                    'blackbody'*     Blackbody illumination (specify at
%                                     minimum the temperature)
%                    'd50'            D50
%                    'd65'            D65
%                    'equal energy'   Equal energy/Uniform/White
%                    'equal photons'  Equal photons
%                    'fluorescent'    Fluorescent
%                    'illuminant c'   Illuminant C
%                    'tungsten'       Tungsten
%                 Some of the optional parameters in the structure are:
%                    name             The name of the illuminant source
%                    luminance        The strength of the illuminant
%                    temperature      The light temperature (required for a
%                                     blackbody illuminant)
%                    spectrum.wave    The wavelength(s)
%    colorSpace - (Optional) The color space for the data. Default is XYZ.
%                 Options for the color space are:
%                    'LAB'            CIELAB
%                    'XYZ'            CIE XYZ
%                    'lRGB'           linear RGB (from sRGB)
%                    'sRGB'           display sRGB
%                    'Stockman'       (not yet implemented)
%                    'SmithPokorny'   (not yet implemented)
%
% Outputs:
%    target     - The macbeth color chart values
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    macbethCompareIdeal, macbethReadReflectance
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    01/31/18  jnm  Formatting

% Examples:
%{
    macbethXYZ = macbethIdealColor('d65', 'xyz');
    macbethXYZ = reshape(macbethXYZ, 4, 6, 3);
    macbethsRGB = xyz2srgb(macbethXYZ);
    vcNewGraphWin;
    image(macbethsRGB);
%}
%{
    macbethXYZ = macbethIdealColor('d65', 'xyz');
    xy = chromaticity(macbethXYZ);
    chromaticityPlot(xy);
%}
%{
    lRGB = macbethIdealColor('d65', 'lrgb');
    lRGB = XW2RGBFormat(lRGB, 4, 6);
    vcNewGraphWin;
    image(lRGB);
%}
%{
    lightParameters.name = 'blackbody';
    lightParameters.temperature = 6000;
    lightParameters.spectrum.wave = 400:10:700;
    lightParameters.luminance = 100;
    lRGB = macbethIdealColor(lightParameters, 'lrgb');
    image(lRGB);
%}
%{
    macbethLAB = macbethIdealColor('tungsten', 'lab');
    plot3(macbethLAB(:, 1), macbethLAB(:, 2), macbethLAB(:, 3), 'ro');
    grid on
    gSeries = 4:4:24;
    hold on
    plot3(macbethLAB(gSeries, 1), macbethLAB(gSeries, 2), ...
        macbethLAB(gSeries, 3), 'kx'); 
    xlabel('L');
    ylabel('a');
    zlabel('b')
%}

if notDefined('illuminant'), illuminant = 'D65'; end
if notDefined('colorSpace'), colorSpace = 'XYZ'; end
if checkfields(illuminant, 'spectrum', 'wave')
    wave = illuminant.spectrum.wave;
else
    wave = (400:10:700)';
end

% For the order of the reflectances, see the comment in macbethChartCreate
% or type: load('macbethChart', 'comment'); comment 
patchList = 1:24;
whitePatch = 4;
macbethChart = macbethReadReflectance(wave, patchList);

% Read illumination.  Could be a string or a structure describing a
% blackbody illuminant.
if ischar(illuminant)
    illEnergy = illuminantRead([], illuminant);
else
    illEnergy = illuminantRead(illuminant);
end

colorSignal = diag(illEnergy) * macbethChart;

switch lower(colorSpace)
    case 'xyz'
        % Sets the max Y value to 100 cd/m2
        macbethXYZ = ieXYZFromEnergy(colorSignal', wave);
        target = 100 * (macbethXYZ / max(macbethXYZ(:, 2)));
    case 'lab'
        macbethXYZ = ieXYZFromEnergy(colorSignal', wave);
        macbethXYZ = 100 * (macbethXYZ / max(macbethXYZ(:, 2)));
        whiteXYZ = macbethXYZ(whitePatch, :);
        target = ieXYZ2LAB(macbethXYZ, whiteXYZ);
    case 'lrgb'
        % Linear RGB values from srgb space; these don't include gamma
        macbethXYZ = macbethIdealColor(illuminant, 'xyz');
        % For xyz2srgb, the XYZ is scaled so that max is around
        % 1. We want the linear RGB to really correspond to the Y values we
        % send in. So, we must scale back.
        macbethXYZ = macbethXYZ / 100;  % Set max to 1 - max Y is 100.
        % vcNewGraphWin; image(xyz2srgb(XW2RGBFormat(macbethXYZ, 4, 6)))
        % Below #ok<ASGLU>
        [idealSRGB, idealLRGB] = xyz2srgb(XW2RGBFormat(macbethXYZ, 1, 24));
        idealLRGB = RGB2XWFormat(idealLRGB);
        
        % Clipping 
        target = ieClip(idealLRGB, 0, 1);
        % vcNewGraphWin; image(XW2RGBFormat(target, 4, 6))
        
    case 'srgb'
        macbethXYZ = macbethIdealColor(illuminant, 'xyz');
        idealSRGB = xyz2srgb(XW2RGBFormat(macbethXYZ, 1, 24));
        target = RGB2XWFormat(idealSRGB);
    case 'stockman'
        disp('Not yet implemented')
    case 'smithpokorny'
        disp('Not yet implemented')
    otherwise
        error('Unknown color space.')
end

end