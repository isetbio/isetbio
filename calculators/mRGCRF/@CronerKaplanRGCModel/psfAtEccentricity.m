function [hEcc, vEcc, thePSFs, thePSFsupportDegs] = psfAtEccentricity(goodSubjects, eccXrange, eccYrange, deltaEcc)

    subtractCentralRefraction = true;
    imposedRefractionErrorDiopters = 1.0
    
    measurementWavelength = 550;
    wavelengthsListToCompute = 550;
    wavefrontSpatialSamples = 501;
    desiredPupilDiamMM = 3;
    
    % For monkey retina
    micronsPerDegree = 300;
        
    for subjIdx = 1:numel(goodSubjects)
        
        subjectIndex = goodSubjects(subjIdx);
        [hEcc, vEcc, zCoeffIndices, zMap, pupilDiamMM] = getTypicalSubjectZcoeffs(subjectIndex, ...
            subtractCentralRefraction, imposedRefractionErrorDiopters, deltaEcc, eccXrange, eccYrange);
       
       for eccYIndex = 1:numel(vEcc)
       for eccXIndex = 1:numel(hEcc)

           fprintf('Computing PSF for subject %d, ecc = %f %f\n', subjectIndex, hEcc(eccXIndex), vEcc(eccYIndex));
            zCoeffs = zeros(1,21);
            zCoeffs(zCoeffIndices+1) = squeeze(zMap(eccYIndex, eccXIndex,:));
            
            % Generate oi at this eccentricity
            theOI = makeCustomOI(zCoeffs, pupilDiamMM, measurementWavelength, ...
                desiredPupilDiamMM, wavelengthsListToCompute, wavefrontSpatialSamples, micronsPerDegree);
            
            % Extract PSF
            [thePSF,thePSFsupportDegs] = extractPSFfromOI(theOI, wavelengthsListToCompute);
            
            % Allocate memory
            if (subjIdx*eccYIndex*eccXIndex == 1)
                thePSFs = zeros(numel(goodSubjects), numel(vEcc), numel(hEcc), size(thePSF,1), size(thePSF,2));
            end
            
            % Save PSF
            thePSFs(subjIdx, eccYIndex, eccXIndex,:,:) = thePSF;
       end
       end
            
    end
%     
end

function [thePSF, thePSFsupportDegs] = extractPSFfromOI(theOI, targetWavelength)
    optics = oiGet(theOI, 'optics');

    wavelengthSupport = opticsGet(optics, 'wave');
    [~,idx] = min(abs(wavelengthSupport-targetWavelength));
    targetWavelength = wavelengthSupport(idx);

    % Get PSF slice at target wavelength
    thePSF = opticsGet(optics,'psf data',targetWavelength);

    % Extract support in arcmin
    psfSupportMicrons = opticsGet(optics,'psf support','um');
    if (isfield(optics, 'micronsPerDegree'))
        micronsPerDegree = optics.micronsPerDegree;
    else
        focalLengthMeters = opticsGet(optics, 'focalLength');
        focalLengthMicrons = focalLengthMeters * 1e6;
        micronsPerDegree = focalLengthMicrons * tand(1);
    end

    xGridDegs = psfSupportMicrons{1}/micronsPerDegree;
    thePSFsupportDegs = xGridDegs(1,:);
end


function theOI = makeCustomOI(zCoeffs, measPupilDiameterMM, measWavelength, ...
    desiredPupilDiamMM, wavelengthsListToCompute, wavefrontSpatialSamples, micronsPerDegree)

    showTranslation = false;
    
    [thePSF, theOTF, xSfCyclesDeg, ySfCyclesDeg, xMinutes, yMinutes, theWVF] = ...
        computePSFandOTF(zCoeffs, ...
             wavelengthsListToCompute, wavefrontSpatialSamples, ...
             measPupilDiameterMM, desiredPupilDiamMM, ...
             measWavelength, showTranslation, 'doNotZeroCenterPSF', true);
        
    for waveIndex = 1:numel(wavelengthsListToCompute)
        theWaveOTF = squeeze(theOTF(:,:,waveIndex));
        theOTF(:,:,waveIndex) = ifftshift(theWaveOTF);
    end
    
    theOI = oiCreate('wvf human', desiredPupilDiamMM,[],[], micronsPerDegree);
    optics = oiGet(theOI,'optics');
    optics = opticsSet(optics, 'otfwave', wavelengthsListToCompute);
    
    % Update optics with new OTF data
    xSfCyclesPerMM = 1000*xSfCyclesDeg / micronsPerDegree;
    ySfCyclesPerMM = 1000*ySfCyclesDeg / micronsPerDegree;
    customOptics = opticsSet(optics,'otf data',theOTF);
    customOptics = opticsSet(customOptics, 'otffx',xSfCyclesPerMM);
    customOptics = opticsSet(customOptics,'otffy',ySfCyclesPerMM);
    
    % Update theOI with custom optics
    theOI = oiSet(theOI,'optics', customOptics);
end

function [hEccQ, vEccQ, zCoeffIndices, zMapQ, pupilDiamMM] = getTypicalSubjectZcoeffs(subjectIndex, subtractCentralRefraction, ...
    refractionErrorDiopters, deltaEcc, ...
    eccXrange, eccYrange)

    minHorizontalEcc = eccXrange(1);
    maxHorizontalEcc = eccXrange(2);
    minVerticalEcc = eccYrange(1);
    maxVerticalEcc = eccYrange(2);
    
     allData = rawDataReadData('zCoefsPolans2015', ...
                    'datatype', 'isetbiomatfileonpath');
     allData = allData.data;
     pupilDiamMM = 4.0;
     
     zCoeffs = squeeze(allData(subjectIndex , :, 3:end));
     zCoeffIndices = 3:size(allData,3);
     
     vEcc = 25:-5:-25;
     hEcc = 40:-1:-40;
     zMap = zeros(numel(vEcc), numel(hEcc),numel(zCoeffIndices));
     pt = 0;
     for vEccIndex = 1:numel(vEcc)
         for hEccIndex = 1:numel(hEcc)
             pt = pt+1;
             zMap(vEccIndex, hEccIndex,:) = zCoeffs(pt,:);
         end
     end
     
     vEccQ = max([minVerticalEcc min(vEcc)]) :deltaEcc: min([max(vEcc) maxVerticalEcc]);
     hEccQ = max([minHorizontalEcc min(hEcc)]) :deltaEcc: min([max(hEcc) maxHorizontalEcc]);
     
     [X,Y] = meshgrid(hEcc, vEcc);
     [XQ,YQ] = meshgrid(hEccQ, vEccQ);
     zMapQ = zeros(numel(vEccQ), numel(hEccQ),numel(zCoeffIndices));
     
     for zIndex = 1:size(zMap,3)
         zz = squeeze(zMap(:,:,zIndex));
         if ((zCoeffIndices(zIndex) == 4) && (subtractCentralRefraction))
             idx = find((X==0) & (Y==0));
             micronsError = wvfDefocusDioptersToMicrons(refractionErrorDiopters, pupilDiamMM);
             zz = zz - zz(idx) + micronsError;
         end
         zz = interp2(X,Y,zz, XQ, YQ);
         zMapQ(:,:,zIndex) = zz;
     end
end