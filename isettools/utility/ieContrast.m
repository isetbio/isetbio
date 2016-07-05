function contrast = ieContrast(sig,varargin)
% Convert a signal to a contrast
%
%    contrast = ieContrast(sig)
%
% The signal mean is subtracted and the signal is scaled so that the entire
% range is 1 unit.  
%
% N.B. If the data are constant, then the return is all zero contrast.
%
% JRG/BW ISETBIO Team, 2016

%% Parse
p = inputParser;
p.addRequired('sig',@isnumeric);

p.parse(sig,varargin{:});
sig = p.Results.sig;

%%
range = max(sig(:)) - min(sig(:));
if range == 0
    warning('Constant data, hence zero contrast.');
    contrast = zeros(size(sig));
else
    mn = mean(sig(:));
    contrast = (sig - mn)/range;
end

end

