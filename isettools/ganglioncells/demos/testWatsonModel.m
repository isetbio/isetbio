function testWatsonModel
    % Unit tests
    RGCmodels.Watson.unitTest.labelValidity();
    RGCmodels.Watson.unitTest.meridianConversions();
    
    % Paper figure plots
    RGCmodels.Watson.plot.figure1();
    RGCmodels.Watson.plot.figure5();
    RGCmodels.Watson.plot.figure9();
    
    % Other demo plots
    spatialSupportUnits = 'mm';   % choose between {'degs', 'mm', 'microns'}
    densityUnits = 'deg^2';         % choose between {'deg^2', and 'mm^2'}
    densityRange = [100 20*1000];   % 
    spatialSupport = struct('maxEcc', 3, 'eccSamples', 50);
    
    
    % Density map for cones
    RGCmodels.Watson.plot.rfDensity2DMaps(100, 'cones', ...
        'spatialSupportUnits', spatialSupportUnits, ...
        'spatialSupport', spatialSupport, ...
        'densityUnits', densityUnits , ...
        'densityRange',densityRange);
    
    % Density map for all ganglion cell types
    RGCmodels.Watson.plot.rfDensity2DMaps(101, 'all ganglion cells', ...
        'spatialSupportUnits', spatialSupportUnits, ...
        'spatialSupport', spatialSupport, ...
        'densityUnits', densityUnits , ...
        'densityRange',densityRange);
    
    % Density map for midget ganglion cells (ON+OFF)
    RGCmodels.Watson.plot.rfDensity2DMaps(102, 'midget ganglion cells', ...
        'spatialSupportUnits', spatialSupportUnits, ...
        'spatialSupport', spatialSupport, ...
        'densityUnits', densityUnits , ...
        'densityRange',densityRange);
    
    fprintf('All unit tests passed\n');
end
