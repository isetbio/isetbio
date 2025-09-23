function [theMaximalExcursionSTFamplitudeSpectrum, ...
    theOptimalOrientation, ...
    theMaximalExcursionSTFphaseSpectrum, ...
    theUnscaledMaximalExcursionSTFamplitudeSpectrum, ...
    theFullSTFamplitudeSpectra, ...
    theFullSTFphaseSpectra, hFigSinusoidalFits] = maximalExcursionSTFfrom2DSTF(orientationDegs, spatialFrequenciesCPD, ...
	spatialPhasesDegs, temporalSupportSeconds, theSTFresponses, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('visualizeFullAndMaximalExcursionSTF', false, @islogical);
    p.addParameter('visualizeSinusoidalFits', false, @islogical);
    p.addParameter('fixedOptimalOrientation', [], @(x)(ischar(x)||isempty(x)||isscalar(x)));
    p.addParameter('axFullSTF',  [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axSTFslice',  [], @(x)(isempty(x)||isa(x, 'handle')));
    p.parse(varargin{:});
    visualizeFullAndMaximalExcursionSTF = p.Results.visualizeFullAndMaximalExcursionSTF;
    visualizeSinusoidalFits = p.Results.visualizeSinusoidalFits;
    fixedOptimalOrientation = p.Results.fixedOptimalOrientation;
    axFullSTF = p.Results.axFullSTF;
    axSTFslice = p.Results.axSTFslice;


    hFigSinusoidalFits = [];

    if (~isempty(fixedOptimalOrientation))

        if (ischar(fixedOptimalOrientation))
            switch (fixedOptimalOrientation)
                case 'OptimalOrientation'
                    for iOri = 1:numel(orientationDegs)
                        [theSTFamplitudeSpectrum, theSTFphaseSpectum] = ...
                            RGCMosaicConstructor.helper.simulateExperiment.stfFromResponseTimeSeries(...
                                spatialFrequenciesCPD, squeeze(theSTFresponses(iOri,:,:)), ...
                                spatialPhasesDegs, temporalSupportSeconds, ...
                                'visualizeSinusoidalFits', visualizeSinusoidalFits);
            
                        theFullSTFamplitudeSpectra(iOri,:) = theSTFamplitudeSpectrum;
                        theFullSTFphaseSpectra(iOri,:) = theSTFphaseSpectum;
                    end % iOri

                     % Normalize to 1 / max amplitude
                    normalizingFactor = 1/max(theFullSTFamplitudeSpectra(:));
                    theFullSTFamplitudeSpectra = theFullSTFamplitudeSpectra * normalizingFactor;
    
                    % Pick the peak response STF as the visual STF for this cell
                    [~,idx] = max(theFullSTFamplitudeSpectra(:));
                    [peakRow, peakCol] = ind2sub(size(theFullSTFamplitudeSpectra), idx);
                    theOptimalOrientation = orientationDegs(peakRow);

                    % Now visualize the fits to the optimal orientation
                    if (visualizeSinusoidalFits)
                        [~,~,hFigSinusoidalFits] = RGCMosaicConstructor.helper.simulateExperiment.stfFromResponseTimeSeries(...
                                spatialFrequenciesCPD, squeeze(theSTFresponses(peakRow,:,:)), ...
                                spatialPhasesDegs, temporalSupportSeconds, ...
                                'visualizeSinusoidalFits', visualizeSinusoidalFits);
                    end

                    theMaximalExcursionSTFamplitudeSpectrum = theFullSTFamplitudeSpectra(peakRow,:);
                    theMaximalExcursionSTFphaseSpectrum = theFullSTFphaseSpectra(peakRow,:);
    
                    theFullSTFamplitudeSpectra = theFullSTFamplitudeSpectra / normalizingFactor;
                    theUnscaledMaximalExcursionSTFamplitudeSpectrum = theMaximalExcursionSTFamplitudeSpectrum / normalizingFactor;
    
                    if (1==2)
                        figure(44);
                        imagesc(spatialFrequenciesCPD, orientationDegs, theFullSTFamplitudeSpectra);
                   
                        hold on;
                        plot(spatialFrequenciesCPD(peakCol), orientationDegs(peakRow), 'kx');
                        axis 'xy'
                        set(gca, 'XScale', 'log', 'XLim', [0.01 100], 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100]);
                        ylabel('orientation (degs)');
                        xlabel('sf (cpd)');
                        pause
                    end
            
                otherwise
                    error('Unknown fixedOptimalOrientation: ''%s''', fixedOptimalOrientation);
            end
        else
            if (isnan(fixedOptimalOrientation))
                % random orientation
                iOri = randi(numel(orientationDegs));
            else
                % Pick the orientation that is closest to the user-specified fixedOptimalOrientation
                [~,iOri] = min(abs(orientationDegs-fixedOptimalOrientation));
                %iOri = find(orientationDegs == fixedOptimalOrientation);
            end
            [theSTFamplitudeSpectrum, theSTFphaseSpectum] = ...
                    RGCMosaicConstructor.helper.simulateExperiment.stfFromResponseTimeSeries(...
                        spatialFrequenciesCPD, squeeze(theSTFresponses(iOri,:,:)), ...
                        spatialPhasesDegs, temporalSupportSeconds, ...
                        'visualizeSinusoidalFits', visualizeSinusoidalFits);
            theFullSTFamplitudeSpectra(1,:) = theSTFamplitudeSpectrum;
            theFullSTFphaseSpectra(1,:) = theSTFphaseSpectum;
                   
            % Normalize to 1 / max amplitude
            normalizingFactor = 1/max(theFullSTFamplitudeSpectra(:));
            theFullSTFamplitudeSpectra = theFullSTFamplitudeSpectra * normalizingFactor;
    
            % Pick the highest extension STF as the visual STF for this cell
            theOptimalOrientation = orientationDegs(iOri);
            theMaximalExcursionSTFamplitudeSpectrum = theFullSTFamplitudeSpectra;
            theMaximalExcursionSTFphaseSpectrum = theFullSTFphaseSpectra;
    
            if (visualizeSinusoidalFits)
                [~,~,hFigSinusoidalFits] = RGCMosaicConstructor.helper.simulateExperiment.stfFromResponseTimeSeries(...
                    spatialFrequenciesCPD, squeeze(theSTFresponses(iOri,:,:)), ...
                    spatialPhasesDegs, temporalSupportSeconds, ...
                    'visualizeSinusoidalFits', visualizeSinusoidalFits);
            end

            theFullSTFamplitudeSpectra = theFullSTFamplitudeSpectra / normalizingFactor;
            theUnscaledMaximalExcursionSTFamplitudeSpectrum = theMaximalExcursionSTFamplitudeSpectrum / normalizingFactor;
    
            for iOri = 1:numel(orientationDegs)
                [theSTFamplitudeSpectrum, theSTFphaseSpectum] = ...
                    RGCMosaicConstructor.helper.simulateExperiment.stfFromResponseTimeSeries(...
                        spatialFrequenciesCPD, squeeze(theSTFresponses(iOri,:,:)), ...
                        spatialPhasesDegs, temporalSupportSeconds, ...
                        'visualizeSinusoidalFits', visualizeSinusoidalFits);
    
                theFullSTFamplitudeSpectra(iOri,:) = theSTFamplitudeSpectrum;
                theFullSTFphaseSpectra(iOri,:) = theSTFphaseSpectum;
            end % iOri
        end
        
    else
    	for iOri = 1:numel(orientationDegs)
            [theSTFamplitudeSpectrum, theSTFphaseSpectum] = ...
                RGCMosaicConstructor.helper.simulateExperiment.stfFromResponseTimeSeries(...
                    spatialFrequenciesCPD, squeeze(theSTFresponses(iOri,:,:)), ...
                    spatialPhasesDegs, temporalSupportSeconds, ...
                    'visualizeSinusoidalFits', visualizeSinusoidalFits&&visualizeFullAndMaximalExcursionSTF);

            theFullSTFamplitudeSpectra(iOri,:) = theSTFamplitudeSpectrum;
            theFullSTFphaseSpectra(iOri,:) = theSTFphaseSpectum;
        end % iOri

        % Normalize to 1 max amplitude
        normalizingFactor = 1/max(theFullSTFamplitudeSpectra(:));
        theFullSTFamplitudeSpectra = theFullSTFamplitudeSpectra * normalizingFactor;

        % Pick the highest extension STF as the visual STF for this cell
        [theOptimalOrientation, theMaximalExcursionSTFamplitudeSpectrum, ...
         theMaximalExcursionSTFphaseSpectrum] = highestExtensionSTF(orientationDegs, spatialFrequenciesCPD, ...
            theFullSTFamplitudeSpectra, theFullSTFphaseSpectra, ...
            fixedOptimalOrientation);

          % Pick the peak response STF as the visual STF for this cell
          iOri = find(orientationDegs == theOptimalOrientation);
          if (visualizeSinusoidalFits)
            [~,~,hFigSinusoidalFits] = RGCMosaicConstructor.helper.simulateExperiment.stfFromResponseTimeSeries(...
                                spatialFrequenciesCPD, squeeze(theSTFresponses(iOri,:,:)), ...
                                spatialPhasesDegs, temporalSupportSeconds, ...
                                'visualizeSinusoidalFits', visualizeSinusoidalFits);
          end

        if (visualizeFullAndMaximalExcursionSTF)

            if (isempty(axFullSTF)) || (isempty(axSTFslice))
                hFig = figure(); clf;
                set(hFig, 'Position', [700 10 490 1170], 'Color', [1 1 1]);

                axFullSTF = subplot(2,1,1);
                axSTFslice = subplot(2,1,2);
            end

            imagesc(axFullSTF,spatialFrequenciesCPD, orientationDegs, log10(theFullSTFamplitudeSpectra));
            axis(axFullSTF, 'square'); axis(axFullSTF, 'xy')
            colormap(axFullSTF, gray);
            set(axFullSTF, 'CLim', [-2 0], 'XTick', spatialFrequenciesCPD, 'YTick', orientationDegs);
            xtickangle(axFullSTF, 0);
            set(axFullSTF, 'FontSize', 16)
            title(axFullSTF, 'FullSTF');
            xlabel(axFullSTF, 'spatial frequency (c/deg)');
            ylabel(axFullSTF, 'orientation (degs)');

            plot(axSTFslice,spatialFrequenciesCPD, theFullSTFamplitudeSpectra, 'k-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5); hold on;
            plot(axSTFslice,spatialFrequenciesCPD, theMaximalExcursionSTFamplitudeSpectrum, 'ro-', 'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 1.0);
            axis(axSTFslice, 'xy')
            set(axSTFslice, 'YLim', [1e-2 1e0], 'YScale', 'log', 'XScale', 'log');
            title(axSTFslice, sprintf('STF slice to be fitted (ORI=%d)', theOptimalOrientation));
            xtickangle(axSTFslice, 0);
            set(axSTFslice, 'FontSize', 16)
            xlabel(axSTFslice,'spatial frequency (c/deg)');
            ylabel(axSTFslice, 'amplitude');
        end

        theFullSTFamplitudeSpectra = theFullSTFamplitudeSpectra / normalizingFactor;
        theUnscaledMaximalExcursionSTFamplitudeSpectrum = theMaximalExcursionSTFamplitudeSpectrum / normalizingFactor;
    end
end

function [theOptimalOrientation, theOptimalSTFMagnitudeSpectrum, theOptimalSTFphaseSpectrum] = highestExtensionSTF(...
	      			orientationsTested, spatialFrequenciesTested, ...
	      			theMeasuredSTFs, theMeasuredSTFphases, ...
                    fixedOptimalOrientation)

    % Normalize to max
    theNormalizedSTFs = theMeasuredSTFs / max(theMeasuredSTFs(:));

    if (~isempty(fixedOptimalOrientation))
        bestOri = find(orientationsTested == fixedOptimalOrientation);

        if (~isempty(bestOri))
            ;
        else
            fprintf('Failed to set optimal orientation as per users request. Specified orientation (%d) was not one of the tested ones\n', fixedOptimalOrientation);
        end
    else
        % Determine the optimal orientation. Just pick the orientation for which the STF
        % extends to the highest SF at half max
        highResSF = linspace(spatialFrequenciesTested(1), spatialFrequenciesTested(end), 100);

        halfMaxSTF = 0.5*max(theNormalizedSTFs(:));
        sfAtHalfMax = zeros(1,numel(orientationsTested));
        for iOri = 1:numel(orientationsTested)
            theSTFatThisOrientation = squeeze(theNormalizedSTFs(iOri,:));
            theSTFatThisOrientation = interp1(spatialFrequenciesTested, theSTFatThisOrientation, highResSF);
            [maxSTF, sfIndexOfMaxSTF] = max(theSTFatThisOrientation);
            sfIndicesToSearch = sfIndexOfMaxSTF:numel(highResSF);
            [~,idx] = min(abs(theSTFatThisOrientation(sfIndicesToSearch)-halfMaxSTF));
            sfAtHalfMax(iOri) = highResSF(sfIndicesToSearch(idx));
        end % iORI
        [~,bestOri] = max(sfAtHalfMax);
    end

    theOptimalOrientation = orientationsTested(bestOri);
    theOptimalSTFMagnitudeSpectrum = theMeasuredSTFs(bestOri,:);
    theOptimalSTFphaseSpectrum = theMeasuredSTFphases(bestOri,:);
end