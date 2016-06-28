function hwr = ieHwrect(signal,effectiveZero)
% Half-wave rectification of the signal data
%
% The effective zero level for the rectification can change, so we send in
% a paramter.  Default is 0.
%
%    hwr = ieHwrect(signal,-50)
%
% JRG/HJ/BW Isetbio Team, 2016

p = inputParser;
p.addRequired('signal',@isnumeric);
p.addOptional('effectiveZero',0,@isnumeric);
p.parse(signal,effectiveZero);

signal = p.Results.signal;
effectiveZero = p.Results.effectiveZero;

%% Do it

hwr = max(signal,effectiveZero);

end
