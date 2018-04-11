function p = identityLine(ax)
% Draw an identity line on the current axis
%
% Syntax:
%   p = identityLine(ax)
%
% Description:
%    Draw an identity line on the current axis
%
% Inputs:
%    ax - The axis in question
%
% Outputs:
%    p  - The line to be drawn
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/xx       (c) Stanford VISTA Team
%    12/15/17  jnm  Formatting
%    01/19/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    plot(1:10, randn(1, 10), 'o')
    identityLine(gca);
%}

if notDefined('ax'), ax = gca; end

% Minimum and maximum of axes
xlim = get(ax, 'xlim');
ylim = get(ax, 'ylim');
mn = min(xlim(1), ylim(1));
mx = max(xlim(2), ylim(2));

% Here's the line
p = line([mn mx], [mn mx], 'color', [.5 .5 .5], 'linestyle', '--');

% Set line properties.  These probably want to come in as an argument
set(p, 'linewidth', 2);
grid on

end