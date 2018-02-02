function t_rankWVFsubjects()
    
    loadAndProcessRawData = false;
    
    if (loadAndProcessRawData)
        % WVF properties
        p.wavelengthsListToCompute = [450 500 550 600];
        p.wavefrontSpatialSamples = 501;
        p.pupilDiamMM = 3.0;

        % Trimmed ranges
        p.psfRange = 5*[-1 1];
        p.otfRange = [0.5 100];

        % Centering params
        p.centeringWavelength = 550;
        p.centeringProperty = 'mass';      % choose between {'max', 'mass'}

        [allPSFs, allOTFs] = loadTrimAndCenterData(p);
    else
        % Load previously processed data
        load('allOTFandPSFdata.mat', 'allPSFs', 'allOTFs', 'p');
    end
    
    
    projectionSpaceDim = 200;     % how many SVD component to use  
    distanceMetric = 'RMS';   % select between {'RMS', 'ABSDIFF', 'MAX'}
    [sortedSubjects, distances] = sortSubjectsBasedOnSVDDistanceToReference(allPSFs, projectionSpaceDim, distanceMetric);
    
    wavelengthsNum = size(allPSFs,3);
    for wIndex = 1:wavelengthsNum
        plotSortedSubjectData(wIndex, sortedSubjects, distances, squeeze(allPSFs(:,:,wIndex)), p.wavelengthsListToCompute(wIndex));
    end
end

function plotSortedSubjectData(figNo, sortedSubjects, distances, allPSFs, wavelength)
    Npixels = sqrt(size(allPSFs,2));
    
    colsNum = 20;
    rowsNum = 10;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', rowsNum, ...
       'colsNum', colsNum, ...
       'heightMargin',   0.01, ...
       'widthMargin',    0.01, ...
       'leftMargin',     0.01, ...
       'rightMargin',    0.00, ...
       'bottomMargin',   0.01, ...
       'topMargin',      0.01);
   
    hFig = figure(figNo); clf;
    set(hFig, 'Name', sprintf('%dnm', wavelength), 'Position', [10 10 2200 1170]);
    subplot('Position', subplotPosVectors(1,1).v);
    imagesc(reshape(squeeze(allPSFs(1,:)), [Npixels Npixels]));
    title('mean z-coeff');
    axis 'square';
    drawnow;
    
    for k = 2:200
        kk = k-1;
        subjectIndex = sortedSubjects(kk);
        row = floor((k-1)/colsNum)+1;
        col = mod(k-1,colsNum)+1;
        subplot('Position', subplotPosVectors(row,col).v);
        imagesc(reshape(squeeze(allPSFs(subjectIndex,:)), [Npixels Npixels]));
        title(sprintf('# %d D=%2.3f', subjectIndex, distances(subjectIndex)));
        axis 'square';
        set(gca, 'XTick', [], 'YTick', []);
        drawnow;
     end
     colormap(jet(512));
     
end

function [sortedSubjects, distances] = sortSubjectsBasedOnSVDDistanceToReference(allPSFs, projectionSpaceDim, distanceMetric)

    % Get dimensions
    subjectsNum = size(allPSFs,1);
    wavelengthsNum = size(allPSFs,3);
    
    % Normalize all PSFs to unit max
    for sIndex = 1:subjectsNum
        for wIndex = 1:wavelengthsNum
            allPSFs(sIndex,:, wIndex) = allPSFs(sIndex,:, wIndex)/max(squeeze(allPSFs(sIndex,:, wIndex)));
        end
    end
    
    % compute distances for each wavelength
    for wIndex = 1:wavelengthsNum
        wavePSFs = squeeze(allPSFs(:,:,wIndex));
        
        [u, s, v] = svd(wavePSFs, 'econ');
        
        %reconstructedPSFs = u * s *v';
        %varianceExplained = (diag(s)).^2;
         
        for subjectIndex = 1:subjectsNum
            inputVector = squeeze(wavePSFs(subjectIndex,:));
            basisVector = squeeze(v(:,1:projectionSpaceDim));
            SVDspaceCoords(subjectIndex, :, wIndex) = inputVector*basisVector/norm(inputVector);
        end
        
        maxSVDcoordAtThisWavelength = max(max(abs(SVDspaceCoords(:, :, wIndex))));
        SVDspaceCoords(:, :, wIndex) = SVDspaceCoords(:, :, wIndex) / maxSVDcoordAtThisWavelength;
    end % wIndex
    
    % Find distances to the meanZcoeff PSF across all wavelengths 
    SVDspaceCoords = reshape(SVDspaceCoords, [subjectsNum projectionSpaceDim*wavelengthsNum]);
    meanZcoeffSVDspaceCoords = SVDspaceCoords(1,:);
    residuals = bsxfun(@minus, SVDspaceCoords, meanZcoeffSVDspaceCoords);
    
    switch distanceMetric
        case 'RMS'
            distances = sqrt(sum(residuals.^2,2));
        case 'ABSDIFF'
            distances = sum(abs(residuals),2);
        case 'MAX'
            distances = max(abs(residuals),[],2);
        otherwise
            error('Unknown distance metric: ''%s''.', distanceMetric)
    end
    
    [~, sortedSubjects] = sort(distances);    
end

function [allPSFs, allOTFs, xMinutes, yMinutes, xSfCyclesDeg, ySfCyclesDeg] = loadTrimAndCenterData(p)

    % Load Thibos data
    [Zcoeffs_SampleMean, Zcoeffs_S, subject_coeffs] = ...
        wvfLoadThibosVirtualEyes(p.pupilDiamMM);
    ZcoeffSubjects = subject_coeffs.bothEyes;
    subjectsNum = size(ZcoeffSubjects,2);
    
    % Compute the PSFs and OTFs for the sample mean Zcoeffs
    [meanZcoeffPSF, meanZcoeffOTF, xSfCyclesDeg, ySfCyclesDeg, xMinutes, yMinutes] = ...
        computePSFandOTF(Zcoeffs_SampleMean, p.wavelengthsListToCompute, ...
        p.wavefrontSpatialSamples, p.pupilDiamMM, p.centeringWavelength, p.centeringProperty);

    % Trim data (so we can compute SVD later on)
    psfColsToKeep = find(abs(xMinutes) < max(p.psfRange));
    psfRowsToKeep = find(abs(yMinutes) < max(p.psfRange));
    xMinutes = xMinutes(psfColsToKeep);
    yMinutes = yMinutes(psfRowsToKeep);
    
    otfColsToKeep = find(abs(xSfCyclesDeg) < max(p.otfRange));
    otfRowsToKeep = find(abs(ySfCyclesDeg) < max(p.otfRange));
    xSfCyclesDeg = xSfCyclesDeg(otfColsToKeep);
    ySfCyclesDeg = ySfCyclesDeg(otfRowsToKeep );
    
    % Preallocate memory for all PSFs and all OTFs
    allPSFs = zeros(1+subjectsNum, numel(psfRowsToKeep)*numel(psfColsToKeep), numel(p.wavelengthsListToCompute));
    allOTFs = zeros(1+subjectsNum, numel(otfRowsToKeep)*numel(otfColsToKeep), numel(p.wavelengthsListToCompute));
    
    % Add first entry, the mean Z coeff PSF
    allPSFs(1, :, :) = reshape(...
        meanZcoeffPSF(psfRowsToKeep, psfColsToKeep, :), ...
        [1 numel(psfRowsToKeep)*numel(psfColsToKeep) numel(p.wavelengthsListToCompute)]);
    
    % Add first entry, the mean Z coeff OTF
    allOTFs(1, :, :) = reshape(...
        meanZcoeffOTF(otfRowsToKeep, otfColsToKeep, :), ...
        [1 numel(otfRowsToKeep)*numel(otfColsToKeep) numel(p.wavelengthsListToCompute)]);
    
    % Add each subject's data
    for subjectIndex = 1:subjectsNum
        zCoeffs = squeeze(ZcoeffSubjects(:,subjectIndex));
        [subjectPSF, subjectOTF, ~, ~, ~, ~] = ...
            computePSFandOTF(zCoeffs, p.wavelengthsListToCompute, ...
            p.wavefrontSpatialSamples, p.pupilDiamMM, p.centeringWavelength, p.centeringProperty);
        
        % Add subject PSF
        allPSFs(1+subjectIndex, :, :) = reshape(...
            squeeze(subjectPSF(psfRowsToKeep, psfColsToKeep, :)), ...
            [1 numel(psfRowsToKeep)*numel(psfColsToKeep) numel(p.wavelengthsListToCompute)]);
    
        % Add subject OTF
        allOTFs(1+subjectIndex, :, :) = reshape(...
            squeeze(subjectOTF(otfRowsToKeep, otfColsToKeep, :)), ...
            [1 numel(otfRowsToKeep)*numel(otfColsToKeep) numel(p.wavelengthsListToCompute)]);
    end % subjectIndex
    
    save('allOTFandPSFdata.mat', 'allPSFs', 'allOTFs', 'p', '-v7.3');
    fprintf('data saved\n');
end

function [PSFs, OTFs, xSfCyclesDeg, ySfCyclesDeg, xMinutes, yMinutes] = ...
    computePSFandOTF(Zcoeffs, wavelengthsListToCompute, wavefrontSpatialSamples, pupilDiamMM, centeringWavelength, centeringProperty)
    %% Compute WVF
    umPerDegree = 300;
    theWVF = makeWVF(wavefrontSpatialSamples, Zcoeffs, wavelengthsListToCompute, ...
            pupilDiamMM, pupilDiamMM, umPerDegree, '');
    
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
                    translationVector, centeringProperty, ...
                    xSfGridCyclesDegGrid,ySfGridCyclesDegGrid, showTranslation);
        translationVector
    else
        translationVector = [0 0];
    end
    
    for wIndex = 1:numel(wavelengthsListToCompute)
        theWaveOTF = wvfGet(theWVF, 'otf', wavelengthsListToCompute(wIndex));
        theWavePSF = wvfGet(theWVF, 'psf', wavelengthsListToCompute(wIndex));
        showTranslation = false;
        [theWaveOTF, ~, ~, ~,~] = otfWithZeroCenteredPSF(...
            theWaveOTF, theWavePSF, ...
            translationVector, centeringProperty,...
            xSfGridCyclesDegGrid, ySfGridCyclesDegGrid, showTranslation);
        
        [xGridMinutes,yGridMinutes, theWavePSF] = ...
            OtfToPsf(xSfGridCyclesDegGrid,ySfGridCyclesDegGrid,fftshift(theWaveOTF));
        if (wIndex == 1)
            OTFs = zeros(size(theWaveOTF,1), size(theWaveOTF,2), numel(wavelengthsListToCompute));
            PSFs = zeros(size(theWavePSF,1), size(theWavePSF,2), numel(wavelengthsListToCompute));
        end
        
        OTFs(:,:,wIndex) = fftshift(abs(theWaveOTF));
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

function [centeredOTF,  translationVector, centeredPSF, xGridMinutes,yGridMinutes] = ...
    otfWithZeroCenteredPSF(OTF, PSF, translationVector,  centeringProperty, xSfGridCyclesDegGrid, ySfGridCyclesDegGrid, showTranslation)

    if (isempty(translationVector))
        if (strcmp(centeringProperty, 'mass'))
            centerOfPSF = computeCenterOfMass(PSF);
        elseif (strcmp(centeringProperty, 'max'))
            centerOfPSF = computeCenterOfMax(PSF);
        else
            error('Unknown center property: ''%s''.', centeringProperty);
        end
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
        translationVector = (centerPosition-centerOfPSF);
    end
    centeredOTF = shiftInFTplane(OTF, translationVector);

    [xGridMinutes,yGridMinutes,centeredPSF] = OtfToPsf(xSfGridCyclesDegGrid,ySfGridCyclesDegGrid,fftshift(centeredOTF));
    
    if (abs(sum(centeredPSF(:))-1) > 1000*eps(1))
        fprintf('centeredPSF min = %1.9f, max = %1.9f, sum = %1.9f\n', min(centeredPSF(:)), max(centeredPSF(:)), sum(centeredPSF(:)));
        error('PSF volume does not equal 1\n');
    end
   
end

function centerOfMass = computeCenterOfMass(PSF)
    [rc,cc] = ndgrid(1:size(PSF,1),1:size(PSF,2));
    Mt = sum(PSF(:));
    centerOfMassY = sum(PSF(:) .* rc(:)) / Mt;
    centerOfMassX = sum(PSF(:) .* cc(:)) / Mt;
    centerOfMass = [centerOfMassX centerOfMassY];
end

function centerOfMax = computeCenterOfMax(PSF)
    [~,idx] = max(PSF(:));
    [i,j] = ind2sub(size(PSF), idx);
    centerOfMax = [j i];
end

function otf = shiftInFTplane(otf, translationVector)
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