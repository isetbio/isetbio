function visualizeArtalOptics()

    subjectIndices = 1:15;
    whichEyes = {'left eye', 'right eye'};
    micronsPerDegree = 281;
    subtractCentralRefraction = false;
    zeroCenterPSF = true;
    flipPSFUpsideDown = ~true;
    horizEcc = -20:1:20;
    pupilDiameterMM = 4;
    
    sv = NicePlot.getSubPlotPosVectors(...
           'colsNum', numel(horizEcc), ...
           'rowsNum', numel(subjectIndices), ...
           'heightMargin',  0.03, ...
           'widthMargin',    0.003, ...
           'leftMargin',     0.01, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.03, ...
           'topMargin',      0.02); 
       
     for subjectGroup = 1:9
        subjectIndices = (1:15) + (subjectGroup-1)*15;
        for eyeIndex = 1:numel(whichEyes)

            hFig = figure(subjectGroup*10 + eyeIndex);
            clf;
            set(hFig, 'Position', [10 10 2829 1304]);
            cMap = brewermap(1025, '*spectral');

            colormap(cMap);
            for k = 1:numel(subjectIndices)
                if (subjectIndices(k) > 130)
                    continue;
                end

                for xEcc = 1:numel(horizEcc)
                     wavelengthsListToCompute = 550;
                     [theOI, thePSF, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength] = ArtalOptics.oiForSubjectAtEccentricity(...
                         subjectIndices(k), whichEyes{eyeIndex}, [horizEcc(xEcc) 0], ...
                         pupilDiameterMM, wavelengthsListToCompute, micronsPerDegree, ...
                         'subtractCentralRefraction', subtractCentralRefraction, ...
                         'wavefrontSpatialSamples', 501, ...
                         'zeroCenterPSF', zeroCenterPSF, ...
                         'flipPSFUpsideDown', flipPSFUpsideDown);
                     subplot('Position', sv(k,xEcc).v);
                     imagesc(psfSupportMinutesX, psfSupportMinutesY, thePSF);
                     axis 'image';
                     set(gca, 'XLim', 10*[-1 1], 'YLim', 10*[-1 1], 'XTick', [], 'YTick', []);
                     if (k == 1)
                            title(sprintf('%2.0f', horizEcc(xEcc)));
                     end
                     if (xEcc == 1)
                         if (strcmp(whichEyes{eyeIndex}, 'left eye'))
                            ylabel(sprintf('subj. %2.0f (LE)',subjectIndices(k)));
                         else
                             ylabel(sprintf('subj. %2.0f (RE)',subjectIndices(k)));
                         end

                     end
                     drawnow;
                end
            end
            
            % Export Figure
            NicePlot.exportFigToPDF(sprintf('group_%d_%s.pdf', subjectGroup,whichEyes{eyeIndex}), hFig, 300);
        end
     end
    
end
