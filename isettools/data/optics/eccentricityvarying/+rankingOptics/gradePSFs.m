function gradePSFs(exportsDir, analyzedOpticsParams)
	xPos = analyzedOpticsParams.eccentricityPosDegs.x;
	yPos = analyzedOpticsParams.eccentricityPosDegs.y;
    switch (analyzedOpticsParams.whichEye)
        case 'right eye'
	       dataFile = sprintf('%sOpticsRE_AnalysisEccPosDegs_%2.1f_%2.1f.mat', ...
                analyzedOpticsParams.zernikeDataBase, xPos, yPos);
           case 'left eye'
           dataFile = sprintf('%sOpticsLE_AnalysisEccPosDegs_%2.1f_%2.1f.mat', ...
                analyzedOpticsParams.zernikeDataBase, xPos, yPos);
    end

	% Allocate memory
	subjectPSFData = cell(1, max(analyzedOpticsParams.subjectIndices));
	
	for iSubject = 1:numel(analyzedOpticsParams.subjectIndices)

		fprintf('Analyzing PSF %d of %d\n', iSubject, numel(analyzedOpticsParams.subjectIndices));
		analyzedOpticsParams.subjectID = analyzedOpticsParams.subjectIndices(iSubject);
        analyzedOpticsParams.subtractCentralRefraction = analyzedOpticsParams.subjectRequiresCentralRefractionCorrection(iSubject);
                
        % Compute PSFimage and cone image
    	[psfImage, coneMosaicImage, NyquistFrequency, psfSupportArcMin, theZCoeffs] = psfAndConeImage([xPos yPos], analyzedOpticsParams);
               

        hfig = figure(10);
        ax = subplot(1,2,1);
        imagesc(ax, psfSupportArcMin, psfSupportArcMin, psfImage);
        axis(ax, 'image'); axis(ax, 'xy');
        set(ax, 'XLim', 5*[-1 1], 'YLim', 5*[-1 1]);
        set(ax, 'FontSize', 16);
        xlabel(ax, 'arc min'); ylabel(ax, 'arc min');
        title(ax, 'PSF');

        ax = subplot(1,2,2);
        imagesc(ax, psfSupportArcMin, psfSupportArcMin, coneMosaicImage);
        axis(ax, 'image'); axis(ax, 'xy');
        set(ax, 'XLim', 5*[-1 1], 'YLim', 5*[-1 1]);
        set(ax, 'FontSize', 16);
        xlabel(ax, 'arc min'); ylabel(ax, 'arc min');
        title(ax, 'cone aperture');

        colormap(gray(1024));
        drawnow;

        % Compute power spectra
        % The computed cutoff frequencies are the -15 dB amplitude corner frequencies. % At -15dB, the MTF power is 3.16%, and its amplitude is 17.78%
        % cornerFrequencyAttenuationDB = -15;
        % cornerFrequencyAttenuation = 10^(cornerFrequencyAttenuationDB/20);

        [psfImageSpectrum, psfImageSpectrumRoatated, coneMosaicImageSpectrum, ...
         spectralSupportCyclesPerDegree, spectralSupportCyclesPerDegreePositive, coneSpectrumSlice, ...
         spectralSupportCyclesPerDegreePositiveX, psfSpectrumSliceX, ...
         spectralSupportCyclesPerDegreePositiveY, psfSpectrumSliceY, ...
         effectivePSFSpectrumSliceX, effectivePSFSpectrumSliceY, ...
         coneCutoffSF, psfXCutoffSF, psfYCutoffSF] = computePSFAndConePowerSpectra(psfImage, coneMosaicImage, psfSupportArcMin);

    	% Data for exporting
        subjectPSFData{iSubject} = struct(...
            'subjectID', analyzedOpticsParams.subjectID, ...
            'zCoeffs', theZCoeffs, ...
            'whichEye', analyzedOpticsParams.whichEye, ...
            'pupilDiamMM', analyzedOpticsParams.pupilDiameterMM, ...
            'psfXCutoffSF', psfXCutoffSF, ...
            'psfYCutoffSF', psfYCutoffSF, ...
            'PSFSpectrumSliceX', effectivePSFSpectrumSliceX, ...
            'PSFSpectrumSliceY', effectivePSFSpectrumSliceY, ...
            'spectralSupportCyclesPerDegreeX', spectralSupportCyclesPerDegreePositiveX, ...
            'spectralSupportCyclesPerDegreeY', spectralSupportCyclesPerDegreePositiveY, ...
            'coneSpectrumSlice', coneSpectrumSlice, ...
            'psfImage', psfImage, ...
            'coneMosaicImage', coneMosaicImage, ...
            'psfSupportArcMin', psfSupportArcMin);
	end % for iSubject

	% Export data
	save(fullfile(exportsDir, dataFile), 'subjectPSFData');
end
