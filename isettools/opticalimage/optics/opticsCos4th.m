function oi = opticsCos4th(oi)
% Compute relative illumination for cos4th model
%
% Syntax:
%   oi = opticsCos4th(oi)
%
% Description:
%    This routine is used for shift-invariant optics, when full ray trace
%    information is unavailable.
%
%    There are examples contained in the code. To access, type 'edit
%    opticsCos4th.m' into the Command Window.
%
% Inputs:
%    oi - Struct. The optical image structure.
%
% Outputs:
%    oi - Struct. The modified optical image structure.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    03/08/18  jnm  Formatting
%    04/07/18  dhb  Fixed example.

% Examples:
%{
    scene = sceneCreate;
    oi = oiCreate('human');
    oi = oiCompute(oi,scene);
    oi = opticsCos4th(oi)
%}

optics = oiGet(oi, 'optics');
method = opticsGet(optics, 'cos4thfunction');
if isempty(method)
    method = 'cos4th';
    oi = oiSet(oi, 'optics cos4thfunction', optics);
end

% Calculating the cos4th scaling factors
%
% We might check whether it exists already and only do this if the cos4th
% slot is empty.
optics = feval(method, optics, oi);
% figure; mesh(optics.cos4th.value)
oi = oiSet(oi, 'optics', optics);

% Applying cos4th scaling.
sFactor = opticsGet(optics,'cos4thData');  % figure(3); mesh(sFactor)
photons = bsxfun(@times, oiGet(oi, 'photons'), sFactor);

% Compress the calculated image and put it back in the structure.
oi = oiSet(oi, 'photons',photons);

end