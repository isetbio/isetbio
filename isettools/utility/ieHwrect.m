function hwr = ieHwrect(signal, effectiveZero)
% Half-wave rectification of the signal data
%
% Syntax:
%   hwr = ieHwrect(signal, effectiveZero)
%
% Description:
%    The effective zero level for rectification can depend on the signal we
%    are analyzing. This routine lets the user  specify the effective zero
%    level as a parameter (Default is 0). Any signal values below
%    effectiveZero  are set to effectiveZero.
%
% Inputs:
%    signal        - The wave to rectify
%    effectiveZero - The effective zero level to apply to the signal
%
% Outputs:
%    hwr           - Half-wave rectification
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/16  JRG/HJ/BW  Isetbio Team, 2016
%    11/22/17  jnm  Formatting
%    01/16/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    wv = -10:10;
    hwr = ieHwrect(wv, 0)
%}
p = inputParser;
p.addRequired('signal', @isnumeric);
p.addOptional('effectiveZero', 0, @isnumeric);
p.parse(signal, effectiveZero);

signal = p.Results.signal;
effectiveZero = p.Results.effectiveZero;

%% Do it
hwr = max(signal, effectiveZero);

end
