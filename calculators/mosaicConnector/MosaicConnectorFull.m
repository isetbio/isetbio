function MosaicConnectorFull

    rootDir = fileparts(which(mfilename()));
    tmpDir = fullfile(rootDir, 'tmpMatFiles');
    exportsDir = fullfile(rootDir, 'exports');
    
    doInitialMosaicCropping = true;                    % phase 1 - crop within circular window
    checkMosaicSeparationAndCropAgain = true;          % phase 2 - check separation and possibly crop within rectangular window
    assignConeTypes = true;                            % phase 3 - assign cone types
    connectConesToRGCcenters = true;                    % phase 4 - connect cones to RGC RF centers
    visualizeConeToRGCcenterConnections = true;         % phase 5
    coVisualizeRFsizeWithDendriticFieldSize = ~true;     % Phase 6
    
    % Configure the phase run parameters
    connector = containers.Map();
    
    inputMosaic = struct('fov',20, 'eccSamples', 256);
    inputMosaic = struct('fov',40, 'eccSamples', 482);
    
    % Phase1: isolate the central (roiRadiusDeg) mosaic
    connector('phase1') = struct( ...
        'run', doInitialMosaicCropping, ...
        'whichEye', 'right', ...                // input mosaic params
        'mosaicFOVDegs', inputMosaic.fov, ...                // input mosaic params
        'eccentricitySamplesNumCones', inputMosaic.eccSamples, ... // input mosaic params
        'eccentricitySamplesNumRGC', inputMosaic.eccSamples, ...   // input mosaic params
        'roiRadiusDeg', Inf, ...                 // processing params
        'outputFile','roiMosaic', ...
        'outputDir', tmpDir ...
    );
        
    % Phase 2: Check that rfs are not less than a treshold (x mean spacing). If any rfs
    % are found with a separation less than that, an error is thrown. Also
    % crop to desired size
    roiCropDegs = struct('xo', 0.0, 'yo', 0.0, 'width', Inf, 'height', 1.5);
    postFix = sprintf('eccX_%2.1f_eccWidth_%2.1f_eccHeight_%2.1f', roiCropDegs.xo, roiCropDegs.width, roiCropDegs.height);
    
    connector('phase2') = struct( ...
        'run', checkMosaicSeparationAndCropAgain, ...
        'inputFile', connector('phase1').outputFile, ...
        'roiRectDegs', roiCropDegs, ...                                 // cropped region to include
        'thresholdFractionForMosaicIncosistencyCorrection', 0.6, ...    // separation threshold for error
        'outputFile', sprintf('%s__CheckedStats', sprintf('%s_%s',connector('phase1').outputFile, postFix)),...
        'outputDir', tmpDir ...
    );

    % Phase 3: Assign cone types to the coneMosaic
    connector('phase3') = struct( ...
        'run', assignConeTypes, ...
        'inputFile', connector('phase2').outputFile, ...
        'tritanopicAreaDiameterMicrons', 2.0*1000*WatsonRGCModel.rhoDegsToMMs(0.3/2), ...        // Tritanopic region size: 0.3 deg diam
        'relativeSconeSpacing', 2.7, ...                        // This results to around 8-9% S-cones
        'LtoMratio', 2.0, ...                                   // Valid range: [0 - Inf]
        'outputFile', sprintf('%s__ConesAssigned', connector('phase2').outputFile),...
        'outputDir', tmpDir ...
    );

    % Phase 4: Connect cones to the mRGC RF centers
    connector('phase4') = struct( ...
        'run', connectConesToRGCcenters, ...
        'inputFile', connector('phase3').outputFile, ...
        'orphanRGCpolicy', 'remove', ...                        // How to deal with RGCs that have no input
        'outputFile', sprintf('%s__ConeConnectionsToRGCcenters', connector('phase3').outputFile),...
        'outputDir', tmpDir ...
    );


    % Phase 5: Visualize connections
    connector('phase5') = struct( ...
        'run', visualizeConeToRGCcenterConnections, ...
        'inputFile', connector('phase4').outputFile, ...
        'zLevels', [0.3 1], ...                                     // contour levels
        'whichLevelsToContour', [1], ...                            // Which level to plot the isoresponse contour
        'displayEllipseInsteadOfContour', ~true, ...                // Ellipse fitted to the sum of flat-top gaussians
        'patchEccDegs', [18 0], ...                                 // Eccenticity of visualized patch
        'patchSizeDegs', [1.5 1.5], ...                             // Size of visualized patch
        'outputFile', sprintf('%s__Visualization', connector('phase4').outputFile),...
        'exportsDir', exportsDir, ...
        'outputDir', tmpDir ...
    );
        
    % Phase 6: Visualize relationship between RF sizes and dendritic size data
    connector('phase6') = struct( ...
        'run', coVisualizeRFsizeWithDendriticFieldSize, ...
        'inputFile', connector('phase4').outputFile, ...
        'zLevels', [0.3 1], ...                                     // contour levels
        'whichLevelsToContour', [1], ...                            // Which level to fit the ellipse to
        'patchEccDegs', [12 0], ...                                 // Eccenticity of visualized patch
        'patchSizeDegs', [2 2], ... 
        'outputFile', sprintf('%s__Visualization2', connector('phase4').outputFile),...
        'exportsDir', exportsDir, ...
        'outputDir', tmpDir ...
    );

        
    if (connector('phase1').run)
        runPhase1(connector('phase1'));
        %return;
    end
       
    if (connector('phase2').run)
        runPhase2(connector('phase2'));
        %return;
    end
    
    if (connector('phase3').run)
        runPhase3(connector('phase3'));
        %return;
    end
    
    if (connector('phase4').run)
        runPhase4(connector('phase4'));
        %return;
    end
    
    if (connector('phase5').run)
        runPhase5(connector('phase5'));
        %return;
    end
    
    if (connector('phase6').run)
        runPhase6(connector('phase6'));
        %return;
    end   
    
end
