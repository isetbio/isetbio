function [theCustomOI, Zcoeffs, theWVF] = oiWithCustomOptics(opticsModel, wavefrontSpatialSamples, calcPupilDiameterMM, umPerDegree, varargin)
    p = inputParser;
    p.addParameter('showTranslation', false, @islogical);
    p.addParameter('centeringWavelength', 550, @isnumeric);
    p.addParameter('wavelengths', [], @isnumeric);
    p.parse(varargin{:});
    
    showTranslation = p.Results.showTranslation;
    
    % Save rng state. 
    currentState = rng;
    
    % Set rng seed to 1 to ensure we always get the same Zcoeffs
    rng(1);

    % Generate default OI
    theOI = oiCreate('wvf human', calcPupilDiameterMM,[],[], umPerDegree);
    optics = oiGet(theOI,'optics');
    
    if (isempty(p.Results.wavelengths))
        wavelengthsListToCompute = opticsGet(optics,'wave');
    else
        wavelengthsListToCompute = p.Results.wavelengths;
        optics = opticsSet(optics, 'otfwave', wavelengthsListToCompute);
    end
    
    if (strfind(opticsModel, 'AOoptics'))
        fprintf('Generating wavefront object for ''%s'' optics \n', opticsModel);
        [Zcoeffs_SampleMean, ~, ~] = wvfLoadThibosVirtualEyes(7.5);
        Zcoeffs = Zcoeffs_SampleMean * 0;
        measPupilDiameterMM = calcPupilDiameterMM;
        % Generate the WVF    
        theWVF = makeWVF(wavefrontSpatialSamples, Zcoeffs, wavelengthsListToCompute, ...
            measPupilDiameterMM, calcPupilDiameterMM, umPerDegree, opticsModel);
        % Compute an optics structure from the WVF
        optics = oiGet(wvf2oi(theWVF),'optics');
        theCustomOI = oiSet(theOI,'optics', optics);   
        return;
    end
    
    % Load appropriate Thibos dataset
    availableThibosMeasPupilSizes = [3 4.5 6 7.5];
    idx = find(availableThibosMeasPupilSizes>= calcPupilDiameterMM);
    if (isempty(idx))
        error('There are no Thibos data for pupil size large enough for %2.2f mm', calcPupilDiameterMM);
    else
        measPupilDiameterMM = availableThibosMeasPupilSizes(idx(1));
        fprintf('Using the %3.2f mm Thibos data set to compute optics for %3.2f mm pupil.\n', measPupilDiameterMM, calcPupilDiameterMM);
    end
    
    [Zcoeffs_SampleMean, Zcoeffs_S, subject_coeffs] = wvfLoadThibosVirtualEyes(measPupilDiameterMM);
    ZcoeffSubjects = subject_coeffs.bothEyes;
       
    switch opticsModel
        case {'WvfHumanMeanOTFmagMeanOTFphase', 'WvfHumanMeanOTFmagZeroOTFphase'}
            subjectsNum = 100;
            Zcoeffs = ieMvnrnd(Zcoeffs_SampleMean, Zcoeffs_S, subjectsNum)'; 

        case 'WvfHuman'
            Zcoeffs = Zcoeffs_SampleMean;
            
        case 'ThibosBestPSFSubject3MMPupil'
            Zcoeffs = ZcoeffSubjects(:,98);

        case 'ThibosDefaultSubject3MMPupil'
            Zcoeffs = ZcoeffSubjects(:,132);
                
        case 'ThibosAverageSubject3MMPupil'
            Zcoeffs = ZcoeffSubjects(:,54);
                
        case 'ThibosDefocusedSubject3MMPupil'
            Zcoeffs = ZcoeffSubjects(:,194);
                
        case 'ThibosVeryDefocusedSubject3MMPupil'
            Zcoeffs = ZcoeffSubjects(:,21);
                
        otherwise
            error('Unknown optics Model: ''%s''.', opticsModel);
    end
    
    subjectsNum = size(Zcoeffs,2);
    for subjectIndex = 1:subjectsNum
        if (subjectsNum>1)
            fprintf('Generating wavefront object for ''%s'' optics (subject %d of %d)\n', opticsModel, subjectIndex,subjectsNum);
        else
            fprintf('Generating wavefront object for ''%s'' optics \n', opticsModel);
        end
        
        [thePSF, theOTF, xSfCyclesDeg, ySfCyclesDeg, xMinutes, yMinutes, theWVF] = computePSFandOTF(...
            squeeze(Zcoeffs(:,subjectIndex)), ...
            wavelengthsListToCompute, wavefrontSpatialSamples, calcPupilDiameterMM, ...
            p.Results.centeringWavelength, showTranslation);
        
        if (subjectIndex == 1)
            theMeanOTF = theOTF;
        else
            theMeanOTF = theMeanOTF + theOTF;
        end
    end % subjectIndex
    
    if (subjectsNum>1)
        meanOTF = theMeanOTF/subjectsNum;
        
        % Reconstruct from mag and phase
        for waveIndex = 1:numel(wavelengthsListToCompute)
            theWaveOTF = squeeze(meanOTF(:,:,waveIndex));

            otfMeanMag = abs(theWaveOTF);
            otfMeanPhase = angle(theWaveOTF);

            phaseToUse = otfMeanPhase;
            if strcmp(opticsModel, 'WvfHumanMeanOTFmagZeroOTFphase')
                phaseToUse = 0*otfMeanPhase;
            end
        
            theWaveOTF = otfMeanMag .* exp(1i*phaseToUse);
            meanOTF(:,:,waveIndex) = ifftshift(theWaveOTF);
        end
    
        theOTF = meanOTF;
    else
        for waveIndex = 1:numel(wavelengthsListToCompute)
            theWaveOTF = squeeze(theOTF(:,:,waveIndex));
            theOTF(:,:,waveIndex) = ifftshift(theWaveOTF);
        end
    end
    
    % Update optics with new OTF data
    xSfCyclesPerMM = 1000*xSfCyclesDeg / umPerDegree;
    ySfCyclesPerMM = 1000*ySfCyclesDeg / umPerDegree;
    customOptics = opticsSet(optics,'otf data',theOTF);
    customOptics = opticsSet(customOptics, 'otffx',xSfCyclesPerMM);
    customOptics = opticsSet(customOptics,'otffy',ySfCyclesPerMM);
    
    % Update customOI with custom optics
    theCustomOI = oiSet(theOI,'optics', customOptics);
        
    % Restore current rng state
    rng(currentState);
end

function [PSFs, OTFs, xSfCyclesDeg, ySfCyclesDeg, xMinutes, yMinutes, theWVF] = computePSFandOTF(Zcoeffs, wavelengthsListToCompute, wavefrontSpatialSamples, targetPupilDiamMM, centeringWavelength, showTranslation)
    %% Compute WVF
    umPerDegree = 300;
    theWVF = makeWVF(wavefrontSpatialSamples, Zcoeffs, wavelengthsListToCompute, ...
            targetPupilDiamMM, targetPupilDiamMM, umPerDegree, '');
    
    xSfCyclesPerRetinalMicron = wvfGet(theWVF, 'otf support', 'um', wavelengthsListToCompute(1));
    xSfCyclesDeg = xSfCyclesPerRetinalMicron * wvfGet(theWVF,'um per degree');
    ySfCyclesDeg = xSfCyclesDeg;
    [xSfGridCyclesDegGrid,ySfGridCyclesDegGrid] = meshgrid(xSfCyclesDeg, ySfCyclesDeg);
    
    if (~isempty(centeringWavelength))
        % Retrieve OTF at the centeringWavelength
        theCenteringOTF = wvfGet(theWVF, 'otf', centeringWavelength);
        theCenteringPSF = wvfGet(theWVF, 'psf', centeringWavelength);
        translationVector = []; showTranslation = false;
        [~, translationVector, ~, ~, ~] = otfWithZeroCenteredPSF(...
                    theCenteringOTF, theCenteringPSF, ...
                    translationVector, xSfGridCyclesDegGrid,ySfGridCyclesDegGrid, ...
                    showTranslation);
    else
        translationVector = [0 0];
    end
    
    for wIndex = 1:numel(wavelengthsListToCompute)
        theWaveOTF = wvfGet(theWVF, 'otf', wavelengthsListToCompute(wIndex));
        theWavePSF = wvfGet(theWVF, 'psf', wavelengthsListToCompute(wIndex));
        [theWaveOTF, ~, ~, ~,~] = ...
            otfWithZeroCenteredPSF(theWaveOTF, theWavePSF, translationVector,  xSfGridCyclesDegGrid, ySfGridCyclesDegGrid, showTranslation);
        
        [xGridMinutes,yGridMinutes, theWavePSF] = ...
            OtfToPsf(xSfGridCyclesDegGrid,ySfGridCyclesDegGrid,fftshift(theWaveOTF));
        if (wIndex == 1)
            OTFs = zeros(size(theWaveOTF,1), size(theWaveOTF,2), numel(wavelengthsListToCompute));
            PSFs = zeros(size(theWavePSF,1), size(theWavePSF,2), numel(wavelengthsListToCompute));
        end
        
        OTFs(:,:,wIndex) = fftshift(theWaveOTF);
        PSFs(:,:,wIndex) = theWavePSF;
    end
    
    xMinutes = xGridMinutes(1,:);
    yMinutes = yGridMinutes(:,1);
end

function theWVF = makeWVF(wavefrontSpatialSamples, zcoeffs, wavelengthsToCompute, measPupilDiameterMM, calcPupilDiameterMM, umPerDegree, name)
    theWVF = wvfCreate(...
                'umPerDegree', umPerDegree, ...
                'calc wavelengths',wavelengthsToCompute,...
                'measuredpupil', measPupilDiameterMM, ...
                'calc pupil size',calcPupilDiameterMM, ...
                'spatialsamples', wavefrontSpatialSamples, ...
                'zcoeffs', zcoeffs,...
                'name', name);
    
    % Now compute the PSF
    theWVF = wvfComputePSF(theWVF);
end

function [centeredOTF,  translationVector, centeredPSF, xGridMinutes,yGridMinutes] = otfWithZeroCenteredPSF(OTF, PSF, translationVector,  xSfGridCyclesDegGrid, ySfGridCyclesDegGrid, showTranslation)
    if (isempty(translationVector))
        % Compute center of mass
        centerOfMass = computeCenterOfMass(PSF, 'useMaxAsCenterOfMass', true);
        centerPosition = floor(size(PSF,1)/2) + 1 * [1 1];
        
        if (showTranslation)
            figure(200);clf
            imagesc(1:size(PSF,2), 1:size(PSF,1), PSF);
            hold on;
            plot(centerPosition(1)*[1 1], [1 size(PSF,1)], 'k-');
            plot([1 size(PSF,2)], centerPosition(2)*[1 1], 'k-');
            plot(centerOfMass(1), centerOfMass(2), 'ro', 'MarkerFaceColor', [1 0.5 0.5]);
            plot([centerPosition(1) centerOfMass(1)], [centerPosition(2)  centerOfMass(2)], 'r-');
            hold off;
            axis 'square'
            colormap(gray(1024));
            drawnow;
        end
        
        % Compute translation vector to bring center of mass at 0,0
        translationVector = (centerPosition-centerOfMass);
    end
    centeredOTF = shiftInFourierPlane(OTF, translationVector);

    [xGridMinutes,yGridMinutes,centeredPSF] = OtfToPsf(xSfGridCyclesDegGrid,ySfGridCyclesDegGrid,fftshift(centeredOTF));
    
    if (abs(sum(centeredPSF(:))-1) > 1000*eps(1))
        fprintf('centeredPSF min = %1.9f, max = %1.9f, sum = %1.9f\n', min(centeredPSF(:)), max(centeredPSF(:)), sum(centeredPSF(:)));
        error('PSF volume does not equal 1\n');
    end
end

function centerOfMass = computeCenterOfMass(PSF, varargin)
    p = inputParser;
    p.addParameter('useMaxAsCenterOfMass', true, @islogical);
    p.parse(varargin{:});
    
    if (p.Results.useMaxAsCenterOfMass)
        [~,idx] = max(PSF(:));
        [i,j] = ind2sub(size(PSF), idx);
        centerOfMass = [j i];
    else
        [rc,cc] = ndgrid(1:size(PSF,1),1:size(PSF,2));
        Mt = sum(PSF(:));
        centerOfMassY = sum(PSF(:) .* rc(:)) / Mt;
        centerOfMassX = sum(PSF(:) .* cc(:)) / Mt;
        centerOfMass = [centerOfMassX centerOfMassY];
    end
end

function otf = shiftInFourierPlane(otf, translationVector)
    [N, M] = size(otf);
    x = [0:floor(N/2) floor(-N/2)+1:-1];
    y = [0:floor(M/2) floor(-M/2)+1:-1];
    x_shift = exp(-1i * 2 * pi * translationVector(2) * x' / N);
    y_shift = exp(-1i * 2 * pi * translationVector(1) * y / M);

    % Force conjugate symmetry. Otherwise this frequency component has no
    % corresponding negative frequency to cancel out its imaginary part.
    if mod(N, 2) == 0
        x_shift(N/2+1) = real(x_shift(N/2+1));
    end 
    if mod(M, 2) == 0
        y_shift(M/2+1) = real(y_shift(M/2+1));
    end
    maxOTF = max(otf(:));
    otf = otf .* (x_shift * y_shift);
    otf = otf / max(abs(otf(:))) * maxOTF;
end

