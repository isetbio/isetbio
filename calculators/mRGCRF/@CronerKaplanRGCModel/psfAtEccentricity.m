function [hEcc, vEcc, thePSFs, thePSFsupportDegs, theOIs] = psfAtEccentricity(goodSubjects, imposedRefractionErrorDiopters, ...
    desiredPupilDiamMM, wavelengthsListToCompute, micronsPerDegree, wavefrontSpatialSamples, eccXrange, eccYrange, deltaEcc, ...
    varargin)

    % Parse input
    p = inputParser;
    p.addParameter('noLCA', false, @islogical);
    p.addParameter('noOptics', false, @islogical);
    p.parse(varargin{:});
    
    % See if we will simulate no longitudinal chromatic aberration
    noLCAflag = p.Results.noLCA;
    noOpticsFlag = p.Results.noOptics;
    
    subtractCentralRefraction = true;
    % Best focus at 550 nm
    measurementWavelength = 550;
    
    computeMicronsPerDegreeAtEachEccentricity = false;
    if (isempty(micronsPerDegree))
        computeMicronsPerDegreeAtEachEccentricity = true;
    end
    
    for subjIdx = 1:numel(goodSubjects)
        
        subjectIndex = goodSubjects(subjIdx);
        [hEcc, vEcc, zCoeffIndices, zMap, pupilDiamMM] = getTypicalSubjectZcoeffs(subjectIndex, ...
            subtractCentralRefraction, imposedRefractionErrorDiopters, deltaEcc, eccXrange, eccYrange);
       
       for eccYIndex = 1:numel(vEcc)
       for eccXIndex = 1:numel(hEcc)

           fprintf('Computing PSF for subject %d at ecc (x,y) = (%2.1f, %2.1f) degs\n', subjectIndex, hEcc(eccXIndex), vEcc(eccYIndex));
           zCoeffs = zeros(1,21);
           zCoeffs(zCoeffIndices+1) = squeeze(zMap(eccYIndex, eccXIndex,:));
            
           % Compute microns per degree at this eccentricity
           if (computeMicronsPerDegreeAtEachEccentricity)
               eccRadius = sqrt(hEcc(eccXIndex)^2+vEcc(eccYIndex)^2);
               micronsPerDegree = WatsonRGCModel.sizeDegsToSizeRetinalMicrons(1.0, eccRadius);
           end
           
           % Generate oi at this eccentricity
           theOI = makeCustomOIFromPolansSubjectZernikeCoefficients(zCoeffs, pupilDiamMM, measurementWavelength, ...
                desiredPupilDiamMM, wavelengthsListToCompute, wavefrontSpatialSamples, micronsPerDegree, ...
                noLCAflag, noOpticsFlag);
            
           % Extract PSF
           for wIndex = 1:numel(wavelengthsListToCompute)
                [thePSF(:,:,wIndex),thePSFsupportDegs] = extractPSFfromOI(theOI, wavelengthsListToCompute(wIndex));
           end
           
           % Allocate memory
           if (subjIdx*eccYIndex*eccXIndex == 1)
                if (numel(wavelengthsListToCompute)==1)
                    thePSFs = zeros(numel(goodSubjects), numel(vEcc), numel(hEcc), size(thePSF,1), size(thePSF,2));
                else
                    thePSFs = zeros(numel(goodSubjects), numel(vEcc), numel(hEcc), size(thePSF,1), size(thePSF,2), size(thePSF,3));
                end
                theOIs =  cell(numel(goodSubjects), numel(vEcc), numel(hEcc));
           end
            
           % Save PSF and the OI
           if (numel(wavelengthsListToCompute)==1)
                thePSFs(subjIdx, eccYIndex, eccXIndex,:,:) = squeeze(thePSF(:,:,1));
           else
               thePSFs(subjIdx, eccYIndex, eccXIndex,:,:,:) = thePSF;
           end
           theOIs{subjIdx, eccYIndex, eccXIndex} = theOI;
       end
       end
            
    end
%     
end

function [thePSF, thePSFsupportDegs] = extractPSFfromOI(theOI, targetWavelength)
    optics = oiGet(theOI, 'optics');

    wavelengthSupport = opticsGet(optics, 'otfwave');
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


function theOI = makeCustomOIFromPolansSubjectZernikeCoefficients(zCoeffs, measPupilDiameterMM, measWavelength, ...
    desiredPupilDiamMM, wavelengthsListToCompute, wavefrontSpatialSamples, micronsPerDegree, noLCAflag, noOpticsFlag)

    showTranslation = false;
    
    if (noOpticsFlag) 
        fprintf(2, 'Generating Polans optics with no Optics\n');
    else
        if (noLCAflag) 
            fprintf(2, 'Generating Polans optics with no LCA\n');
        else
            fprintf(2, 'Generating standard Polans optics\n');
        end
    end
    
    if (noOpticsFlag)
        zCoeffs = zCoeffs * 0;
    end
    
    [thePSF, theOTF, xSfCyclesDeg, ySfCyclesDeg, xMinutes, yMinutes, theWVF] = ...
        computePSFandOTF(zCoeffs, ...
             wavelengthsListToCompute, wavefrontSpatialSamples, ...
             measPupilDiameterMM, desiredPupilDiamMM, ...
             measWavelength, showTranslation, ...
             'doNotZeroCenterPSF', true, ...
             'micronsPerDegree', micronsPerDegree);
        
    for waveIndex = 1:numel(wavelengthsListToCompute)
        theWaveOTF = squeeze(theOTF(:,:,waveIndex));
        theOTF(:,:,waveIndex) = ifftshift(theWaveOTF);
    end
    
    if (noLCAflag)
        inFocusWindex = find(wavelengthsListToCompute == measWavelength);
        assert(~isempty(inFocusWindex), 'In focus wavelength is not in the list of wavelenghts to compute');
        for waveIndex = 1:numel(wavelengthsListToCompute)
             theOTF(:,:,waveIndex) = theOTF(:,:,inFocusWindex);
        end
    end
    
    theOI = oiCreate('wvf human', desiredPupilDiamMM,[],wavelengthsListToCompute, micronsPerDegree);
    optics = oiGet(theOI,'optics');
    optics = opticsSet(optics, 'wave', wavelengthsListToCompute);
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
    
    % Measurement pupil size (according to Polans paper)
    pupilDiamMM = 4.0;
     
    % Reported Z-coeffs are Z3-Z20
    zCoeffs = squeeze(allData(subjectIndex , :, 3:end));
    % Z3-Z20
    zCoeffIndices = 3:size(allData,3);
     
    % Measured eccentricities
    vEcc = 25:-5:-25;
    hEcc = 40:-1:-40;
    zMap = zeros(numel(vEcc), numel(hEcc),numel(zCoeffIndices));
   
    spatialPointIndex = 0;
    for vEccIndex = 1:numel(vEcc)
         for hEccIndex = 1:numel(hEcc)
             spatialPointIndex = spatialPointIndex+1;
             zMap(vEccIndex, hEccIndex,:) = zCoeffs(spatialPointIndex,:);
         end
    end
    
    % Interpolate measured zMap according to desired eccRange and deltaEcc
    vEccQ = max([minVerticalEcc min(vEcc)]) :deltaEcc: min([max(vEcc) maxVerticalEcc]);
    hEccQ = max([minHorizontalEcc min(hEcc)]) :deltaEcc: min([max(hEcc) maxHorizontalEcc]);
     
    [X,Y] = meshgrid(hEcc, vEcc);
    [XQ,YQ] = meshgrid(hEccQ, vEccQ);
    zMapQ = zeros(numel(vEccQ), numel(hEccQ),numel(zCoeffIndices));
     
    for zIndex = 1:size(zMap,3)
         zz = squeeze(zMap(:,:,zIndex));
         % The 4-th z-coeff is defocus. Subtract central defocus from all
         % spatial positions
         if ((zCoeffIndices(zIndex) == 4) && (subtractCentralRefraction))
             idx = find((X==0) & (Y==0));
             zz = zz - zz(idx);
         end
         
         % Add a uniform refraction defocus error at all positions instread
         % of the measured defocus. This assumes that we are focused on
         % each cell, but with a defocus error.
         if ((zCoeffIndices(zIndex) == 4) && (refractionErrorDiopters ~= 0))
            micronsError = wvfDefocusDioptersToMicrons(refractionErrorDiopters, pupilDiamMM);
            zz = zz*0 + micronsError;
         end
         
         zz = interp2(X,Y,zz, XQ, YQ);
         zMapQ(:,:,zIndex) = zz;
     end
end