function hwr = ieHwrect(signal, effectiveZero)
% Half-wave rectification of the signal data
%
% Syntax:
%   hwr = ieHwrect(signal, effectiveZero)
%
% Description:
%    The effective zero level for the rectification can change, so we send
%    in a parameter. Default is 0. This converts any values below
%    effectiveZero to being effectiveZero.
%
% Inputs:
%    signal        - The wave to rectify
%    effectiveZero - The effective zero level to apply to the signal
%
% Outputs:
%    hwr           - Half-wave rectification
%

% History:
%    xx/xx/16  JRG/HJ/BW  Isetbio Team, 2016
%    11/22/17  jnm  Formatting

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
