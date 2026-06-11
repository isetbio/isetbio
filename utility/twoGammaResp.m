function [resp, defs] = twoGammaResp(t, param)
% Create temporal response composed of the difference of 2 gamma functions
%
% Syntax:
%   [resp, defs] = twoGammaResp(timePoints, param)
%
% Description:
%    Create temporal response composed of the difference of 2 gamma
%    functions. This is still being used for defining the coupling and
%    feedback (cpTR and fbTR). See the new function 'gammaPDF' which might
%    be a replacement for this.
%
%    A gamma function is calculate by this formula:
%       1 / (gamma(gamma)*b)  *  (tP/b).^5  .*  exp(-tP/b);
%
%    The final output is a normalized wieghted sum of the 2 gammas.
%       resp = twoGammaResp(100) -- will return a response of length 100
%
%	 The following will return the same response, and the function defaults
%	 will be contained in defs
%       [resp defs] = twoGammaResp(100)
%
%    The following will return the response, with a different linear
%    scaling based on defs.normF
%       defs.normF = 800;
%       resp = twoGammaResp(100, defs)
%
% Inputs:
%    t     - The time points, in milliseconds.
%    param - A structure containing key properties, such as:
%       b     - The location of the gamma peaks. Default is nT ./ [32 16]).
%               i.e 32% and 16% of the total time length.
%       gamma - The parameter to choose the gamma base. Default is [5 5].
%               This is used in multiplicative coefficients only, and the
%               same for the two function.
%       f     - The factor for each gamma function. Default is [1 .5].
%               (resp = f(1) * gamma1 - f(2) * gamma2).
%       normF - The factor used to enlarge the response. Default is 0.0025.
%               (resp = resp / norm(resp) * normF)
%
% Outputs:
%    resp  - The temporal gamma response.
%    defs  - The function defaults, as designated by param.
%
% Notes:
%    * [Note: XXX - The sum of the returned response is not necessarily 1.
%      In the layer call to create the impulse response of the center and
%      surround we scale the return to be unit so that they preserve the
%      mean. Perhaps that should happen here?]
%    * [Note: JW/BW - Should have units of milliseconds we think. We need
%      to make sure the units are right somehow.]
%    * [Note: EC - What is the use of the gamma field in the params
%      variable? (see inputs section).]
%

% History:
%    xx/xx/09       (c) 2009 Stanford Vista Group
%    12/13/17  jnm  Formatting
%    01/19/18  jnm  Formatting update to match Wiki. Fix outputs section.


%% Set defaults - time units are probably ms at this point
nT = length(t);
defs.b = nT ./ [32 16];     % Peak parameters for the gamma functions
defs.f = [1 .5];            % Factor for each gamma function
defs.gamma = [5 5];             % Gamma factors

%% Check inputs and overwrite defaults

if ~notDefined('param')
    names = fieldnames(param); 
else
    names = [];
end

for ii = 1:length(names)
    try
        defs.(names{ii}) = param.(names{ii});
    catch
        error('Unrecognized input: %s', names{ii});
    end
end

%% Extract values for calculations
b = defs.b;
f = defs.f;
g = defs.gamma;

%% Calculate the 2-gamma functions

% Create the two individual gamma-functions - units are milliseconds. 
% Confirm.
resp1 = 1 / (gamma(g(1)) * b(1)) * (t / b(1)) .^ 5 .* exp(-t / b(1));
resp2 = 1 / (gamma(g(2)) * b(2)) * (t / b(2)) .^ 5 .* exp(-t / b(2));

% Create the full function
resp = f(1) *  resp1 - f(2) * resp2;

% Normalize response to coincide with stimulus intensity
resp = (resp ./ norm(resp));

end