function model = osAdaptSteadyState(obj, bgR, varargin)
% Steady-state background current calculated from the background rates.
%
% Syntax:
%   obj = osAdaptSteadyState(obj, bgR, sz)
%
% Description:
%    This function is called internally from @osBioPhys/osCompute.m
%
%    Steady-state background current calculated from the background rates.
%
%    Examples are contained in the code. To access, type 'edit
%    osAdaptSteadyState.m' into the Command Window.
%
% Inputs:
%    obj      - osBioPhys object
%    bgR      - Vector(/matrix) of background isomerization rates (R */sec)
%    varargin - (Optional) Additional information/parameters.
%
% Output:
%    model    - osBioPhys model parameters needed for computation of
%               adapted current over time series input. The model is
%               attached to the input object and returned as a parameter.
%
% Optional key/value pairs:
%    None.
%
% References:
%    http://isetbio.org/cones/adaptation%20model%20-%20rieke.pdf
%    https://github.com/isetbio/isetbio/wiki/Cone-Adaptation
%
% Notes:
% * [Note: JNM - The example doesn't work, despite my edits. Current
%      error is "Undefined function 'osAdaptSteadyState' for input
%      arguments of type 'osLinear'."]
% * [TODO: TO DISCUSS - What's this about? - from Calculation below.]
% * Programming note: Notice that the computation is a search over a
%      bounded variable. The upper and lower bounds are huge for current
%      and thus adequate to always find a solution (we think).

% History:
%    xx/xx/14  HJ   ISETBIO Team, 2014
%    08/xx/16  JRG  ISETBIO Team, updated 8/2016
%    02/14/18  jnm  Formatting
%    04/05/18  NPC  Fixed broken example

% Examples:
%{
    os1 = osCreate('biophys');
    bgR = 10;
    pRate = .1;
    osAdaptSteadyState(os1, bgR, ...
        [size(pRate, 1) size(pRate, 2)]);
%}

%% Parameters
if notDefined('bgR'), error('Background isomerization rate required.'); end
if isscalar(bgR) && ~isempty(varargin), bgR = bgR * ones(varargin{1}); end

%% Calculation
%  In most cases, the input bgR matrix will only contain only a small
%  number of unique values (the background is almost flat) and we can just
%  loop over them
% [TODO: TO DISCUSS - What's this about?]
% Also, what are the units of bgR?  R*/sec?  Or R*/timeSample?

sz = size(bgR);
[bgR, ~, recover_index] = unique(bgR(:));

q = 2 * obj.model.beta * obj.model.cdark ...
    / (obj.model.k * obj.model.gdark ^ obj.model.h);
smax = obj.model.eta / obj.model.phi * obj.model.gdark ...
    * (1 + (obj.model.cdark / obj.model.kGc) ^ obj.model.n);

model = obj.model;

% Init steady state current for background
bgCur = zeros(length(bgR));
for ii = 1 : length(bgR)
    v = bgR(ii);
    bgCur(ii) = fminbnd(@(x) abs(x - (obj.model.k * obj.model.beta * ...
        obj.model.cdark) * (smax * obj.model.phi) ^ obj.model.h / (v / ...
        obj.model.sigma + obj.model.eta) ^ obj.model.h / ...
        (obj.model.beta * obj.model.cdark + q * x) / (1 + (q * x / ...
        obj.model.beta / obj.model.kGc) ^ obj.model.n) ^ obj.model.h), ...
        0, 1000);
end

% recover background isomerization rate and adapted steady state current to
% the input matrix format
bgR = reshape(bgR(recover_index), sz);
bgCur = reshape(bgCur(recover_index), sz);

% Compute additional initial values
model.opsin = bgR / obj.model.sigma;
model.PDE = (model.opsin + obj.model.eta) / obj.model.phi;
model.Ca = bgCur * q / obj.model.beta;
model.Ca_slow = model.Ca;
model.st = smax ./ (1 + (model.Ca / obj.model.kGc) .^ obj.model.n);
model.cGMP = model.st * obj.model.phi ./ (model.opsin + obj.model.eta);

model.bgCur = bgCur;

end