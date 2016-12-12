function oiFrame = frameAtIndex(obj, index)
% Compute the oi at desired index
%
% We should store the photons as doubles separately.  I asked not to
% duplicate, but in fact the double(oi.data.photons) call is very slow, and
% we shouldn't keep doing it with the get in this routine.
%
% NP/BW ISETBIO Team, 2016

% Extract the fixed and modulated photons 
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


if (strcmp(obj.composition, 'add'))
    retinalPhotons = fixedPhotons + obj.modulationFunction(index)*modulatedPhotons;
else
    retinalPhotons = fixedPhotons*(1-obj.modulationFunction(index)) + obj.modulationFunction(index)*modulatedPhotons;
end

if (~isnan(obj.modulationRegion.radiusInMicrons))
    % modulate a subregion only
    pos = oiGet(obj.oiModulated, 'spatial support', 'microns');
    ecc = sqrt(sum(pos.^2, 3));
    mask = ecc < obj.modulationRegion.radiusInMicrons;
    
    for k = 1:size(retinalPhotons,3)
        fullFrame = retinalPhotons(:,:, k);
        background = fixedPhotons;
        background = background(:,:, k);
        fullFrame(mask == 0) = background(mask == 0);
        retinalPhotons(:,:, k) = fullFrame;
    end
end

% Return the composite oi
oiFrame = obj.oiFixed;
oiFrame = oiSet(oiFrame, 'photons', retinalPhotons);
end
