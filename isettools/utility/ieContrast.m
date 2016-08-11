function contrast = ieContrast(sigIn,varargin)
% Convert a signal to a contrast
%
%    contrast = ieContrast(sig)
%
% The signal mean is subtracted and the signal is scaled so that the entire
% range is 1 unit.  
%
% N.B. If the data are constant, then the return is all zero contrast.
%
% See also irComputeLinearSTSeparable.m
% 
% JRG/BW ISETBIO Team, 2016

%% Parse
p = inputParser;
p.addRequired('sigIn',@isnumeric);
p.parse(sigIn,varargin{:});
sigIn = p.Results.sigIn;

sizeSig = size(sigIn);
sigSS = sigIn;

% % Only compute denominator for contrast using steady state signal
% sigSS = sigIn(:,:,end-round(.5*sizeSig(3)):end);

%%
range = max(sigSS(:)) - min(sigSS(:));
sig = sigIn;
if range == 0
    warning('Constant data, hence zero contrast.');
    contrast = zeros(size(sig));
else
    mn = mean(sig(:));
    contrast = (sig - mn)/range;
end

end

