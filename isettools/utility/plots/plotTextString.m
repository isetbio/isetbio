function t = plotTextString(str, position, delta, fontSize)
% Add a text string to a 2D graph in a specific position.
%
% Syntax:
%   t = plotTextString(str, [position], [delta=0.2], [fontSize=12])
%
% Description:
%    Place a text string on a 2D graph in one of several canonical
%    positions. The background of the text is set to white to make the text
%    visible even if the grid is turned on
%
%    This routine could be generalized to 3D, but it has not yet been.
%
%    Possible positions are:  'ul', 'ur', 'll', 'lr' for upper left, upper
%    right, lower left, and lower right.
%
% Inputs:
%    str      - The text to display on the graph
%    position - (Optional) The position (relative to the graph at which to
%               display the text. Default is 'ur' (upper right). The
%               options are as follows:
%           'ul' - Upper Left
%           'ur' - Upper Right (Default)
%           'll' - Lower Left
%           'lr' - Lower Right
%    delta    - (Optional) The amount to offset text away from the edge by.
%               Default is 0.2
%    fontSize - (Optional) The text's font size. Default is 12.
%
% Outputs:
%    t        - The text struct
%
% Notes:
%    * [Note: XXX - TODO: Positions aren't right. Fix.]
%    * [Note: XXX - TODO: We need to account for the scale type when
%      setting these positions.  Not done properly yet.  Also, it would be
%      better to account for the string length, too.  At the very least, we
%      could count the number of letters to set the value of delta.]
%    * [Note: JNM - TODO: We should add a check for if delta is a single
%      argument vs a pair of arguments. Ex. if length(delta) == 1,
%      delta = [delta delta]; end]
%

% History:
%    xx/xx/06       Copyright Imageval Consulting, LLC 2006
%    12/11/17  jnm  Formatting

% Example:
%{
	txt = 'Hello World';
	t = plotTextString(txt, 'ul');
%}

if notDefined('position'), position = 'ur'; end
if notDefined('delta'), delta = [0.2, 0.2]; end
if notDefined('fontSize'), fontSize = 12; end

xlim = get(gca, 'xlim');
ylim = get(gca, 'ylim');

xscale = get(gca, 'xscale');
yscale = get(gca, 'yscale');

switch lower(position)
    case 'ul'
        x = xlim(1) + (xlim(2) - xlim(1)) * delta(1); 
        y = ylim(2) - (ylim(2) - ylim(1)) * delta(2);
    case 'll'
        x = xlim(1) + (xlim(2) - xlim(1)) * delta(1); 
        y = ylim(1) + (ylim(2) - ylim(1)) * delta(2);
    case 'ur'
        x = xlim(2) - (xlim(2) - xlim(1)) * delta(1); 
        y = ylim(2) - (ylim(2) - ylim(1)) * delta(2);
    case 'lr'
        x = xlim(2) - (xlim(2) - xlim(1)) * delta(1); 
        y = ylim(1) + (ylim(2) - ylim(1)) * delta(2);
    otherwise
        error('Unknown position');
end

% Display
t = text(x, y, str);
set(t, 'Background', 'w', 'Fontsize', fontSize);

return;
