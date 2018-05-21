function str = displayDescription(d)
% Text description of the display properties, shown in display window
%
% Syntax:
%   displayDescription(display)
%
% Description:
%    Text description of the display properties, which is shown in the
%    display window.
%
% Inputs:
%    d   - Struct. An ISETBIO display structure.
%
% Outputs:
%    str - String. A string containing the desired information.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/14  HJ   Created May, 2014
%    05/17/18  jnm  Formatting

% Examples:
%{
    d = displayCreate('LCD-Apple');
    str = displayDescription(d)
%}

if notDefined('d'), d = []; end
global vcSESSION

if isempty(d)
    str = 'No display structure';
else
    str = sprintf('Name:\t%s\n', displayGet(d, 'name'));

    wave = displayGet(d, 'wave');
    spacing = displayGet(d, 'binwidth');
    str = addText(str, sprintf('Wave:\t%d:%d:%d nm\n', ...
        min(wave(:)), spacing, max(wave(:))));

    str = addText(str, sprintf('Num. primaries:\t%d\n', ...
        displayGet(d, 'nprimaries')));
    str = addText(str, sprintf('Color bit depth:\t%d\n', ...
        displayGet(d, 'bits')));
    I = displayGet(d, 'main image');
    % Original width/height display (remove if comfortable with changes)
    % str = addText(str, sprintf('Image width: \t%d\nHeight: \t%d', ...
    %    size(I, 2), size(I, 1)));
    str = addText(str, sprintf('Image Width: \t%d\n', size(I, 2)));
    str = addText(str, sprintf('Image Height: \t%d', size(I, 1)));

end
