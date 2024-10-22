function populationPSFvariationWithEccentricity(ZernikeDataBase, whichEye, visualizePSFs)
%
% ZernikeDataBase = 'Polans2015'; 
% whichEye = 'right eye'
% visualizePSFs = false;
% rankingOptics.populationPSFvariationWithEccentricity(ZernikeDataBase, whichEye, visualizePSFs)
%

	assert(ismember(ZernikeDataBase, {'Artal2012', 'Polans2015'}), 'Zenike data must be set to either ''Artal2012'', or ''Polans2015''.');
	assert(ismember(whichEye, {'right eye', 'left eye'}), 'eye must be set to either ''right eye'', or ''left eye''.');

	switch (ZernikeDataBase)
		case 'Artal2012'
			rankedSubjectIDs = ArtalOptics.constants.subjectRanking(whichEye);
			for iSubj = 1:numel(rankedSubjectIDs)
    			subjectRequiresCentralRefractionCorrection(iSubj) = ...
    				ArtalOptics.constants.subjectRequiresCentralRefractionCorrection(whichEye, rankedSubjectIDs(iSubj));
    		end

    		eccXdegs = ArtalOptics.constants.measurementHorizontalEccentricities;
        	eccYdegs = ArtalOptics.constants.measurementVerticalEccentricities; % = 0;

		case 'Polans2015'
			rankedSubjectIDs = PolansOptics.constants.subjectRanking;
			for iSubj = 1:numel(rankedSubjectIDs)
    			subjectRequiresCentralRefractionCorrection(iSubj) = ...
    				PolansOptics.constants.subjectRequiresCentralRefractionCorrection(rankedSubjectIDs(iSubj));
    		end

    		eccXdegs = PolansOptics.constants.measurementHorizontalEccentricities; 
        	eccYdegs = PolansOptics.constants.measurementVerticalEccentricities;
	end

	% Just the 5 best subjects
	rankingsAnalyzed = 1:5; %10; % 1:5;
	rankedSubjectIDs = rankedSubjectIDs(rankingsAnalyzed);
     
	% Range of eccentricities to examine
	temporalEccDegs = [-36:2:-8 -7:1:0];
	eccXdegs = cat(2, temporalEccDegs, fliplr(-temporalEccDegs(1:end-1)));

	eccXYdegs(:,1) = eccXdegs;
	eccXYdegs(:,2) = eccXYdegs(:,1) * 0;

	heightMargin = 0.1;
    bottomMargin =  0.05;
    topMargin = 0.04;
	sv = NicePlot.getSubPlotPosVectors(...
           'colsNum', size(eccXYdegs,1), ...
           'rowsNum', numel(rankedSubjectIDs), ...
           'heightMargin',  heightMargin, ...
           'widthMargin',    0.01, ...
           'leftMargin',     0.02, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   bottomMargin, ...
           'topMargin',      topMargin); 

	zLevels = 0.1:0.2:1.0;
	xLim = 6*[-1 1];
	yLim = 6*[-1 1];
	xTicks = -10:2:10;
	yTicks = -10:2:10;

	

	if (visualizePSFs)
		hFig = figure(200); clf;
	    set(hFig, 'Position', [10 10 1900 1150], 'Color', [1 1 1], 'Name', 'PSF (550)');
	    cLUT = 1-gray(1024);
	    colormap(cLUT);

	    hFigFits = figure(300); clf;
	    set(hFigFits, 'Position', [10 10 1900 1150], 'Color', [1 1 1], 'Name', 'Gaussian fit to PSF (550)');
	    cLUT = 1-gray(1024);
	    colormap(cLUT)
	end

    coneCutoffSF = nan(1, size(eccXYdegs,1));
    psfXCutoffSF = nan(numel(rankedSubjectIDs), size(eccXYdegs,1));
    psfYCutoffSF = nan(numel(rankedSubjectIDs), size(eccXYdegs,1));
    cutoffSF = nan(numel(rankedSubjectIDs), size(eccXYdegs,1));
    minorAxisCutoffSF = nan(numel(rankedSubjectIDs), size(eccXYdegs,1));
	majorAxisCutoffSF = nan(numel(rankedSubjectIDs), size(eccXYdegs,1));


    psfMinorSigma = nan(numel(rankedSubjectIDs), size(eccXYdegs,1));
    psfMajorSigma = nan(numel(rankedSubjectIDs), size(eccXYdegs,1));
    psfSigma = nan(numel(rankedSubjectIDs), size(eccXYdegs,1));
    coneSigma = nan(1, size(eccXYdegs,1));


	for iEcc = 1:size(eccXYdegs,1)
		eccentricityDegsPosition = eccXYdegs(iEcc,:);
		if ( ...
			(eccentricityDegsPosition(1)>= 13) && ...
			(eccentricityDegsPosition(1)<= 18) ...
		    )
			continue;
		end

		for iSubj = 1:numel(rankedSubjectIDs)
			analyzedOpticsParams = struct(...
		        'zernikeDataBase',  ZernikeDataBase, ...
		        'subjectID', rankedSubjectIDs(iSubj), ...
		        'subtractCentralRefraction', subjectRequiresCentralRefractionCorrection(iSubj), ...
		        'zeroCenterPSF', true, ...
		        'flipPSFUpsideDown', true, ...
		        'pupilDiameterMM', 3, ...
		        'whichEye', whichEye);

			% Compute PSF at 550 nm and cone image
	    	[psf550, coneMosaicImage, ~, psfSupportArcMin] = psfAndConeImage(eccXYdegs(iEcc,:), analyzedOpticsParams);

	    	if (iSubj == 1)
	    		% Fit Gaussian to cone aperture
	    		ellipsoidParams = fitGaussianEllipsoid(coneMosaicImage, psfSupportArcMin);
	    		coneSigma(iEcc) = sqrt(ellipsoidParams.xSigma * ellipsoidParams.ySigma);
	    	end

	    	% Fit Gaussian to PSF (not always a good model)
	    	ellipsoidParams = fitGaussianEllipsoid(psf550, psfSupportArcMin);
	    	psfSigma(iSubj,iEcc) = sqrt(ellipsoidParams.xSigma * ellipsoidParams.ySigma);
			psfMinorSigma(iSubj,iEcc) = min([ellipsoidParams.xSigma  ellipsoidParams.ySigma]);
			psfMajorSigma(iSubj,iEcc) = max([ellipsoidParams.xSigma  ellipsoidParams.ySigma]);

	    	% Compute power spectra
	        [psfImageSpectrum, psfImageSpectrumRoatated, coneMosaicImageSpectrum, ...
	         spectralSupportCyclesPerDegree, spectralSupportCyclesPerDegreePositive, coneSpectrumSlice, ...
	         spectralSupportCyclesPerDegreePositiveX, psfSpectrumSliceX, ...
	         spectralSupportCyclesPerDegreePositiveY, psfSpectrumSliceY, ...
	         effectivePSFSpectrumSliceX, effectivePSFSpectrumSliceY, ...
	         coneCutoffSF(iEcc), psfXCutoffSF(iSubj,iEcc), psfYCutoffSF(iSubj,iEcc)] = computePSFAndConePowerSpectra(psf550, coneMosaicImage, psfSupportArcMin);

	         minorAxisCutoffSF(iSubj,iEcc) = min([psfXCutoffSF(iSubj,iEcc) psfYCutoffSF(iSubj,iEcc)]);
	         majorAxisCutoffSF(iSubj,iEcc) = max([psfXCutoffSF(iSubj,iEcc) psfYCutoffSF(iSubj,iEcc)]);
	         cutoffSF(iSubj,iEcc) = sqrt(psfXCutoffSF(iSubj,iEcc) * psfYCutoffSF(iSubj,iEcc)); 

	        if (visualizePSFs)
	        	figure(hFig);
	    		ax = subplot('Position', sv(iSubj, iEcc).v);
		    	contourf(ax, psfSupportArcMin, psfSupportArcMin, psf550/max(psf550(:)), zLevels);
		        hold(ax, 'on');
		        plot(ax, psfSupportArcMin, psfSupportArcMin*0, 'r-', 'LineWidth', 1.0);
		        plot(ax, psfSupportArcMin*0, psfSupportArcMin, 'r-', 'LineWidth', 1.0);
		        hold(ax, 'off');
		        colorbar(ax)
		        axis(ax, 'image'); axis(ax, 'xy'); grid(ax, 'on');
		        set(ax, 'CLim', [0 1], 'ZLim', [0 1], 'XLim', xLim, 'YLim', yLim, 'XTick', xTicks, 'YTick', yTicks);
		        set(ax, 'FontSize', 16);
		        if (iSubj == 1)
		        	title(ax, sprintf('%2.0f, %2.0f (degs)', eccXYdegs(iEcc,1), eccXYdegs(iEcc,2)));
		        end

		        if (iSubj == numel(rankedSubjectIDs))
		        	xlabel(ax, 'arc min'); 
		        else
		        	set(ax, 'XTickLabel', {});
		        end

		        if (iEcc == 1)
		        	ylabel(ax, 'arc min');
		        else
					set(ax, 'YTickLabel', {});
		        end

		        drawnow;

		        figure(hFigFits);
		        ax = subplot('Position', sv(iSubj, iEcc).v);
		    	contourf(ax, psfSupportArcMin, psfSupportArcMin, ellipsoidParams.fittedPSF/max(ellipsoidParams.fittedPSF(:)), zLevels);
		        hold(ax, 'on');
		        plot(ax, psfSupportArcMin, psfSupportArcMin*0, 'r-', 'LineWidth', 1.0);
		        plot(ax, psfSupportArcMin*0, psfSupportArcMin, 'r-', 'LineWidth', 1.0);
		        hold(ax, 'off');
		        colorbar(ax)
		        axis(ax, 'image'); axis(ax, 'xy'); grid(ax, 'on');
		        set(ax, 'CLim', [0 1], 'ZLim', [0 1], 'XLim', xLim, 'YLim', yLim, 'XTick', xTicks, 'YTick', yTicks);
		        set(ax, 'FontSize', 16);
		        if (iSubj == 1)
		        	title(ax, sprintf('%2.0f, %2.0f (degs)', eccXYdegs(iEcc,1), eccXYdegs(iEcc,2)));
		        end

		        if (iSubj == numel(rankedSubjectIDs))
		        	xlabel(ax, 'arc min'); 
		        else
		        	set(ax, 'XTickLabel', {});
		        end

		        if (iEcc == 1)
		        	ylabel(ax, 'arc min');
		        else
					set(ax, 'YTickLabel', {});
		        end

		        drawnow;

	    	end

    	end  % for iSubj
    end % for iEcc


    mRGCcenterColor = [1 0.6 0.0];
    mRGCsurroundColor = [1 0.5 0.4];

    hFig = figure(2); clf;
    set(hFig, 'Position', [10 10 1000 1200], 'Color', [1 1 1]);
    ax = subplot('Position', [0.08 0.05 0.91 0.94]);
    coneResolutionHandle = plot(ax, squeeze(eccXYdegs(:,1)), coneCutoffSF, 'o-', 'Color', [0 0.7 0], 'LineWidth', 1.0, ...
    	'MarkerSize', 14, 'MarkerFaceColor', [0.5 1.0 0.5], 'MarkerEdgeColor', [0 0.6 0]);
    hold(ax, 'on');
    psfResolutionHandle = plot(ax, squeeze(eccXYdegs(:,1)), mean(cutoffSF, 1, 'omitnan'), 'k-', 'LineWidth', 1.5)
	psfMinorAxisResolutionHandle = plot(ax, squeeze(eccXYdegs(:,1)), mean(minorAxisCutoffSF, 1, 'omitnan'), 'r-', 'LineWidth', 1.0)
    psfMajorAxisResolutionHandle = plot(ax, squeeze(eccXYdegs(:,1)), mean(majorAxisCutoffSF, 1, 'omitnan'), 'b-', 'LineWidth', 1.0)
   
    legend(ax, [coneResolutionHandle psfResolutionHandle(1) psfMinorAxisResolutionHandle psfMajorAxisResolutionHandle], ...
    	{'cone mosaic', sprintf('PSF (subj. ranking: %d - %d)', rankingsAnalyzed(1), rankingsAnalyzed(end)), 'PSF (low res)', 'PSF (high res)'});
    set(ax, 'XLim', [-36 36], 'XTick', [-35:5:35]);
    set(ax, 'YLim', [1 200], 'YTick', [1 3 10 30 100 200], 'YScale', 'log');
    set(ax, 'FontSize', 20);
    grid(ax, 'on');
    xtickangle(ax, 0);
    ylabel(ax, 'resolution (c/deg)');
    xlabel(ax, 'eccentricity (degs)');

    hFig = figure(3); clf;
    set(hFig, 'Position', [10 10 1000 1200], 'Color', [1 1 1]);
    ax = subplot('Position', [0.08 0.05 0.91 0.94]);

    coneSigmaHandle = plot(ax, squeeze(eccXYdegs(:,1)), coneSigma, 'o-', 'Color', [0 0.7 0], 'LineWidth', 1.0, ...
    	'MarkerSize', 12, 'MarkerFaceColor', [0.5 1.0 0.5], 'MarkerEdgeColor', [0 0.6 0]);
    hold(ax, 'on');
    psfSigmaHandle = plot(ax, squeeze(eccXYdegs(:,1)), mean(psfSigma, 1), 'ko-', 'LineWidth', 1.5, ...
    	'MarkerSize', 12, 'MarkerFaceColor', [0.75 0.75 0.75], 'MarkerEdgeColor', [0.2 0.2 0.2]);
   
    psfSigmaPercentiles = prctile(psfSigma, [5 95],1);
    plot(ax, squeeze(eccXYdegs(:,1)), psfSigmaPercentiles(1,:), 'k:', 'LineWidth', 1.5);
    plot(ax, squeeze(eccXYdegs(:,1)), psfSigmaPercentiles(2,:), 'k:', 'LineWidth', 1.5);

	psfMinorSigmaHandle = plot(ax, squeeze(eccXYdegs(:,1)), mean(psfMinorSigma, 1), 'r--', 'LineWidth', 1.5);
    psfMajorSigmaHandle = plot(ax, squeeze(eccXYdegs(:,1)), mean(psfMajorSigma, 1), 'b--', 'LineWidth', 1.5);
   
    % Bring in Croner&Kaplan mRGC data
    % The mRGC centers
    [eccDegs, RcDegs] = RGCmodels.CronerKaplan.digitizedData.parvoCenterRadiusAgainstEccentricity;
    % Temporal meridian
    CK95centerSigmaPlotHandle = plot(ax, -eccDegs, RcDegs/sqrt(2.0)*60, 'mo', 'LineWidth', 1.5, 'MarkerSize', 12, ...
    	'MarkerFaceColor', mRGCsurroundColor, 'MarkerEdgeColor', 0.7*mRGCsurroundColor);
    
    % Nasal meridian (exclude optic disk locations)
    idx = find(eccDegs <= 12 | eccDegs >= 20);
    eccDegs = eccDegs(idx);
    RcDegs = RcDegs(idx);
    plot(ax, eccDegs, RcDegs/sqrt(2.0)*60, 'mo', 'LineWidth', 1.5, 'MarkerSize', 12, ...
    	'MarkerFaceColor', mRGCsurroundColor, 'MarkerEdgeColor', 0.7*mRGCsurroundColor);

    % The mRGC surrounds
    % Temporal meridian
    [eccDegs, RsDegs] = RGCmodels.CronerKaplan.digitizedData.parvoSurroundRadiusAgainstEccentricity;
    CK95surroundSigmaPlotHandle = plot(ax, -eccDegs, RsDegs/sqrt(2.0)*60, 'ms', 'LineWidth', 1.5, 'MarkerSize', 16, ...
    	'MarkerFaceColor', mRGCsurroundColor, 'MarkerEdgeColor', 0.7*mRGCsurroundColor);
    % Nasal meridian (exclude optic disk locations)
    idx = find(eccDegs <= 12 | eccDegs >= 20);
    eccDegs = eccDegs(idx);
    RsDegs = RsDegs(idx);
    plot(ax, eccDegs, RsDegs/sqrt(2.0)*60, 'ms', 'LineWidth', 1.5, 'MarkerSize', 16, ...
    	'MarkerFaceColor', mRGCsurroundColor, 'MarkerEdgeColor', 0.7*mRGCsurroundColor);

    legend(ax, [coneSigmaHandle psfSigmaHandle(1) psfMinorSigmaHandle(1) psfMajorSigmaHandle(1) CK95centerSigmaPlotHandle CK95surroundSigmaPlotHandle], ...
    	{'cone mosaic', sprintf('PSF (subj. ranking: %d - %d)', rankingsAnalyzed(1), rankingsAnalyzed(end)), 'PSF (minor axis)', 'PSF (major axis)', 'CK95 mRGC centers' 'CK95 mRGC surrounds'});
    grid(ax, 'on');
    set(ax, 'XLim', [-36 36], 'XTick', [-35:5:35]);
    set(ax, 'YScale', 'log');
    set(ax, 'YLim', [0.1 120], 'YTick', [0.1 0.3 1 3 6 10 15 30 60 120]);
    set(ax, 'FontSize', 20);
    
    xtickangle(ax, 0);
    ylabel(ax, 'sigma (arc min)');
    xlabel(ax, 'eccentricity (degs)');

end

function ellipsoidParams = fitGaussianEllipsoid(psf550, psfSupportArcMin)

	[maxPSF, idx] = max(psf550(:));
	[rowMax, colMax] = ind2sub(size(psf550), idx);
	psf550 = psf550 / maxPSF;

	maxXo = psfSupportArcMin(colMax);
	maxYo = psfSupportArcMin(colMax);
	maxSigmaX = max(psfSupportArcMin);
	maxSigmaY = maxSigmaX;


	%               [gain,   xo,                    xSigma,                 yo,                    ySigma,       orientation]
	initialParams = [1       maxXo                  0.1*maxSigmaX           maxYo                  0.1*maxSigmaY     0.0];
    lowerBounds   = [0.8     min(psfSupportArcMin)  0                       min(psfSupportArcMin)  0               -pi/2];
    upperBounds   = [1.2     max(psfSupportArcMin)  maxSigmaX               max(psfSupportArcMin)  maxSigmaY   		pi/2];


    [X,Y] = meshgrid(psfSupportArcMin, psfSupportArcMin);
	spatialSupport = zeros(size(X,1),size(Y,2),2);
	spatialSupport(:,:,1) = X;
	spatialSupport(:,:,2) = Y;

	options = optimoptions(@lsqcurvefit,'Algorithm','levenberg-marquardt','MaxIterations',1500, 'UseParallel', true);
    [bestFitParams,resnorm,residual,exitflag] = ...
    	lsqcurvefit(@GaussianEllipsoid, initialParams, spatialSupport, psf550, ...
        lowerBounds, upperBounds, options);

    if (exitflag == 1)
    	fprintf('----> Fit failed !!!\n');
    end

    ellipsoidParams.gain = bestFitParams(1);
	ellipsoidParams.x0 = bestFitParams(2);
	ellipsoidParams.xSigma = bestFitParams(3);
	ellipsoidParams.y0 = bestFitParams(4);
	ellipsoidParams.ySigma = bestFitParams(5);
	ellipsoidParams.rotationDegs = bestFitParams(6)/pi*180;

	ellipsoidParams.fittedPSF = GaussianEllipsoid(bestFitParams, spatialSupport);


end

function F = GaussianEllipsoid(params, xdata)	
	gain = params(1);
	x0 = params(2);    
 	y0 = params(4);
	xSigma = params(3); 
	ySigma = params(5);
	theRotation = params(6);

	cosPhi = cos(theRotation);
	sinPhi = sin(theRotation);

	theRotationMatrix = [...
		cosPhi  -sinPhi; ...
		sinPhi   cosPhi];

	X = xdata(:,:,1); Y = xdata(:,:,2);
	n = size(X,1); m = size(X,2);

	XYrot = [X(:) Y(:)] * theRotationMatrix;
	xy0rot = [x0 y0] * theRotationMatrix;

	X = XYrot(:,1) - xy0rot(1);
	Y = XYrot(:,2) - xy0rot(2);

	F = gain * exp(-0.5*(X/xSigma).^2) .* exp(-0.5*(Y/ySigma).^2);
	F = reshape(F, [n m]);
end

