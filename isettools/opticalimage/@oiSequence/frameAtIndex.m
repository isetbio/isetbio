% Compute the oi at desired index
function oiFrame = frameAtIndex(obj, index)

    if (strcmp(obj.composition, 'add'))
        retinalPhotons = obj.fixedPhotons + obj.modulationFunction(index)*obj.modulatedPhotons;
    else
        retinalPhotons = obj.fixedPhotons*(1-obj.modulationFunction(index)) + obj.modulationFunction(index)*obj.modulatedPhotons;
    end

    if (~isnan(obj.modulationRegion.radiusInMicrons))
        % modulate a subregion only
        pos = oiGet(obj.oiModulated, 'spatial support', 'microns');
        ecc = sqrt(sum(pos.^2, 3));
        mask = ecc < obj.modulationRegion.radiusInMicrons;

        for k = 1:size(retinalPhotons,3)
            fullFrame = retinalPhotons(:,:, k);
            background = obj.fixedPhotons;
            background = background(:,:, k);
            fullFrame(mask == 0) = background(mask == 0);
            retinalPhotons(:,:, k) = fullFrame;
        end
    end

    % Return the composite oi
    oiFrame = obj.oiFixed;
    oiFrame = oiSet(oiFrame, 'photons', retinalPhotons);
end
