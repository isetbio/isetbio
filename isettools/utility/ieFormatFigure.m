function ieFormatFigure(fig, fontName, fontSize, figSize, border)
% Formats a figure for presentations or papers.
%
% Syntax:
%   ieFormatFigure([fig], [fontName], [fontSize], [figSize], [border])
%
% Description:
%    The font style, font size, figure size and border size can be adjusted.
%    All parameters are optional. Using this makes it easier to save the
%    figure in an appropriate format for Adobe Illustrator. 
%
% Inputs:
%    This function has no required inputs.
%
%    fig      - (Optional)Figure handle.
%               Default 0
%    fontName - (Optional) Name of the font as a string.
%               Default 'Helvetica'
%    fontSize - Size of the font (points) as a vector.
%               Default [18 14]
%               Format: [axes_labels tick_labels]
%    figSize  - Size of the figure [width height] (inches) as a vector.
%               Default [6 6]
%    border   - Space around the figure (inches) as a vector.
%               Default [0.75 0.35]
%               Format: [left bottom right top] or
%                       [left/right  botom/top]
%
% Outputs:
%    None.
%
% Notes:
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    11/22/17  jnm  Formatting

% Examples:
%{
    ieFormatFigure();
    ieFormatFigure(0,'Arial');
%}

% Format the input parameters.
% Load default paramemters.
if notDefined('fig'), fig = 0; end
if notDefined('fontname'), fontName = 'Helvetica'; end
if notDefined('fontsize'), fontSize = [18 14]; end
if notDefined('figsize'), figSize = [6.5 6.5]; end
if notDefined('border'), border = [1 0.5]; end

% Check the fontsize.
if (length(fontSize) == 1), fontSize = [fontSize fontSize]; end

% Check the figure size.
if (length(figSize) == 1), figSize = [figSize figSize]; end

% Check the border.
if (length(border) == 2)
    border = [border(1) border(1) border(2) border(2)];
end

% Get the figure axes.
if (~notDefined('fig') ~= 1 || fig == 0), fig = gcf; end
axs = get(fig, 'CurrentAxes');

% Get the current figure position.
if (strcmp(get(fig, 'Units'), 'inches') == 1)
    pos = get(fig, 'Position'); 
else
    pos = [0 0 0 0];
end

% Set the figure properties.
set(fig, 'Units', 'inches');
set(fig, 'Position', [pos(1:2) figSize]);
set(fig, 'PaperPosition', [4.25 - figSize(1) / 2, 5.5 - figSize(2) / 2, ...
    figSize]);
set(fig, 'Color', [1 1 1]);

% Set the axes properties.
set(get(axs, 'Title'), 'FontName', fontName, 'FontSize', fontSize(1));
set(get(axs, 'XLabel'), 'FontName', fontName, 'FontSize', fontSize(1));
set(get(axs, 'YLabel'), 'FontName', fontName, 'FontSize', fontSize(1));
set(get(axs, 'ZLabel'), 'FontName', fontName, 'FontSize', fontSize(1));
set(axs, 'FontName', fontName, 'FontSize', fontSize(2));

set(axs, 'Units', 'inches');
set(axs, 'Position', [border(1:2), (figSize - border(1:2) - border(3:4))]);

end
