function theOI = wvf2oiSpecial(theWVF, umPerDegree, pupilDiameterMM)
    % Generate oi from the wvf
    theOI = wvf2oi(theWVF);
    
    % Adjust the OI's fNumber and focalLength to be consistent with the
    % micronsPerDegree and pupilDiameter of the WVF
    optics = oiGet(theOI, 'optics');
    focalLengthMM = (umPerDegree * 1e-3) / (2 * tand(0.5));
    focalLengthMeters = focalLengthMM * 1e-3;

    pupilRadiusMeters = (pupilDiameterMM / 2) * 1e-3;
    pupilDiameterMeters = 2 * pupilRadiusMeters;
    optics = opticsSet(optics, 'fnumber', focalLengthMeters / pupilDiameterMeters);
    optics = opticsSet(optics, 'focalLength', focalLengthMeters);
    theOI = oiSet(theOI, 'optics', optics);
end