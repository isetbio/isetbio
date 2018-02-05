function t = plotTextString(str, position, delta, fontSize)
% Add a text string to a 2D graph in a specific position.
%
% Syntax:
%   t = plotTextString(str, [position], [delta], [fontSize])
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
%    t        - The text string
%
% Optional key/value pairs:
%    None.
%
% Notes:
%

% History:
%    xx/xx/06       Copyright Imageval Consulting, LLC 2006
%    12/11/17  jnm  Formatting
%    12/26/17   BW  Added examples and responded to JNM
%    01/24/18  jnm  Formatting update to match Wiki.

% Example:
%{
	txt = 'Hello World';
    vcNewGraphWin; plot(1:10,1:10,'o');
	t = plotTextString(txt, 'ul');
	t = plotTextString(txt, 'ul',[.1 .1])
	t = plotTextString(txt, 'ul',[.1 .1],18)

	t = plotTextString(txt, 'lr',[.1 .1],18)
	t = plotTextString(txt, 'll',[.1 .1],10)

%}

%%
if notDefined('position'), position = 'ur'; end
if notDefined('delta'), delta = [0.2, 0.2]; end
if notDefined('fontSize'), fontSize = 12; end

xlim = get(gca, 'xlim');
ylim = get(gca, 'ylim');

if length(delta) == 1, delta(2) = delta(1); end

%%
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

%% Display
t = text(x, y, str);
set(t, 'Background', 'w', 'Fontsize', fontSize);

end
