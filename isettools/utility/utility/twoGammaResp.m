function [resp, defs] = twoGammaResp(t, param)
% Create a temporal response composed of the difference of 2 gamma functions
%
% It is still be used for defining the coupling and feedback (cpTR and
% fbTR).  See the new function 'gammaPDF' which might be a replacement for
% this.
%
%  [resp defs] = twoGammaResp(timePoints, param)
%
% twoGammaResp creates a response function, based on the difference between
% two gamma functions. The structure param defines key properties.
%
% t - Should have units of milliseconds we think. We need to make sure the
% units are right somehow (JW/BW) 
% param - 
%  : field      : function                          : default 
%   - b         - the location of the gamma peaks   - nT ./ [32 16] (i.e. at 32% and 16% of the total time length)
%   - gamma     - parameter to choose gamma base    - [5 5]; (used in multiplicative coeff only, and the same for the two function... what is the use of this? -EC)
%   - f         - factor for each gamma function    - [1 .5]; (resp = f(1)*gamma1 - f(2)*gamma2)
%   - normF     - factor to enlarge the response    - .0025; (resp = resp/norm(resp) * normF)
%
% A gamma function is calculate by this formula:
%   1 / (gamma(gamma)*b)  *  (tP/b).^5  .*  exp(-tP/b);
%
% The final output is a normalized wieghted sum of the 2 gammas.
%
% resp = twoGammaResp(100)
%   will return a response of length 100
%
% [resp defs] = twoGammaResp(100)
%   will return the same response, and the function defaults in defs
%
% defs.normF = 800;
% resp = twoGammaResp(100, defs)
%   will return the response, but now with a different linear scaling
%
% NOTE:  The sum of the returned response is not necessarily 1.  In the
% layer call to create the impulse response of the center and surround we
% scale the return to be unit so that they preserve the mean.  Perhaps that
% should happen here?
%
% (c) 2009 Stanford Vista Group


%% Set defaults - time units are probably ms at this point
nT          = length(t);
defs.b      = nT ./ [32 16];     % Peak parameters for the gamma functions
defs.f      = [1 .5];            % Factor for each gamma function
defs.gamma  = [5 5];             % Gamma factors

%% Check inputs and overwrite defaults

if ~notDefined('param'), names = fieldnames(param); 
else                     names = []; 
end

for ii = 1 : length(names)
    try
        defs.(names{ii}) = param.(names{ii});
    catch
        error('Unrecognized input: %s', names{ii});
    end
end

%% Extract values for calculations
b       = defs.b;
f       = defs.f;
g       = defs.gamma;

%% Calculate the 2-gamma functions

% Create the two individual gamma-functions - units are milliseconds. 
% Confirm.
resp1 = 1 / (gamma(g(1))*b(1)) * (t/b(1)).^5 .* exp(-t/b(1));
resp2 = 1 / (gamma(g(2))*b(2)) * (t/b(2)).^5 .* exp(-t/b(2));

% Create the full function
resp = f(1) *  resp1 - f(2) * resp2;

% Normalize response to coincide with stimulus intensity
resp = (resp ./ norm(resp));

end