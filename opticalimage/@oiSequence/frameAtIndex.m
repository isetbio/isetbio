function oiFrame = frameAtIndex(obj, index)
% Compute the oi at desired index
%
% Syntax:
%   oiFrame = frameAtIndex(obj, index)
%
% Description:
%    Compute the oi at the desired index.
%
% Inputs:
%    obj     - Object. The oi sequence object.
%    index   - Numeric. The index at which to calculate the oi.
%
% Outputs:
%    oiFrame - Matrix. The calculated oi frame.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/16  NP/BW  ISETBIO Team, 2016
%    03/27/18  jnm    Formatting
%    07/02/19  JNM    Formatting update

% Extract the fixed and modulated photons.
% We do it this way because oiGet
if ~isempty(obj.photonsFixed)
    fixedPhotons = obj.photonsFixed;
else
    fixedPhotons = oiGet(obj.oiFixed, 'photons');
    obj.photonsFixed = fixedPhotons;
end

if ~isempty(obj.photonsModulated)
    modulatedPhotons = obj.photonsModulated;
else
    modulatedPhotons = oiGet(obj.oiModulated, 'photons');
    obj.photonsModulated = modulatedPhotons;
end

if strcmp(obj.composition, 'add')
    retinalPhotons = fixedPhotons + ...
        obj.modulationFunction(index) * modulatedPhotons;
elseif strcmp(obj.composition, 'blend')
    retinalPhotons = fixedPhotons * (1 - obj.modulationFunction(index)) ...
        + obj.modulationFunction(index) * modulatedPhotons;
elseif strcmp(obj.composition, 'xor')
    if obj.modulationFunction(index) == 0
        retinalPhotons = fixedPhotons;
    else
        retinalPhotons = obj.modulationFunction(index) * modulatedPhotons;
    end
else
    error('Unknown oiSequence composition: ''%s''.', obj.composition);
end

if ~isnan(obj.modulationRegion.radiusInMicrons)
    % modulate a subregion only
    pos = oiGet(obj.oiModulated, 'spatial support', 'microns');
    ecc = sqrt(sum(pos .^ 2, 3));
    mask = ecc < obj.modulationRegion.radiusInMicrons;

    for k = 1:size(retinalPhotons, 3)
        fullFrame = retinalPhotons(:, :, k);
        background = fixedPhotons;
        background = background(:, :, k);
        fullFrame(mask == 0) = background(mask == 0);
        retinalPhotons(:, :, k) = fullFrame;
    end
end

% Return the composite oi
oiFrame = obj.oiFixed;
oiFrame = oiSet(oiFrame, 'photons', retinalPhotons);
end
