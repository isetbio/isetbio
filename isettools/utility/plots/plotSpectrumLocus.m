function [p, l] = plotSpectrumLocus(fig)
% Draw the outline of spectrum locus on the chromacity diagram
%
% Syntax:
%   [p, l] = plotSpectrumLocus([fig])
%
% Description:
%    Draw the outline of the spectrum locus on the chromaticity
%    diagram.  It is a white background with a grid.
%
% Inputs:
%    fig - The figure on which to draw
%
% Outputs:
%    p   - The plot
%    l   - The line
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    chromaticityPlot
%

% History:
%    xx/xx/03       Copyright Imageval 2003
%    12/11/17  jnm  Formatting
%    01/24/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    plotSpectrumLocus;
%}

if notDefined('fig')
    vcNewGraphWin; 
else
    figure(fig);
end

wave = 370:730;
XYZ = ieReadSpectra('XYZ', wave);

% Here are the shifts in the chromaticity of the display
% as the display intensity shifts
spectrumLocus = chromaticity(XYZ);

% These are the (x, y) points of the spectral lines
p = plot(spectrumLocus(:, 1), spectrumLocus(:, 2), '-');
hold on;

% Add a line to close up the outer rim of the spectrum locus curve
l = line([spectrumLocus(1, 1), spectrumLocus(end, 1)], ...
    [spectrumLocus(1, 2), spectrumLocus(end, 2)]);
hold on;
axis equal;
grid on

return;