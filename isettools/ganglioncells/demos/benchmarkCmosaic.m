function benchmarkCmosaic

    % Mosaic options
    varyMPdensity = true;
    varyApertureAndOSlength = true;
    
    doIt('old', varyMPdensity, varyApertureAndOSlength);
    doIt('new', varyMPdensity, varyApertureAndOSlength);
end

function doIt(env, varyMPdensity, varyApertureAndOSlength)
    
    integrationTime = 1000/1000;
    
    % Set data path
    dataPath = fullfile(isetRootPath, 'ganglioncells/demos');
    
    % Compute the oi to the test scene. Also get microns per degree from
    % the oi so we can set it to the cone mosaics.
    [theOI, micronsPerDegreeFromOptics] = computeTheOI(dataPath);

    % For testing - uniform field oi
    %photons = oiGet(theOI, 'photons');
    %photons = bsxfun(@plus, photons*0, mean(mean(photons,1),2));
    %theOI = oiSet(theOI, 'photons', photons);

   
    simulationFOVdegs = 0.5;
    mosaicFOVdegs = simulationFOVdegs*[1 0.5];
    mosaicFOVdegs = [1 1];
    
    if (strcmp(env, 'old'))
        regenerateOldMosaic = true;
        if (regenerateOldMosaic)
            
            theOldConeMosaic = coneMosaicHex(9, ...         % hex lattice sampling factor
                    'fovDegs', mosaicFOVdegs, ...        % match mosaic width to stimulus size
                    'eccBasedConeDensity', true, ...
                    'maxGridAdjustmentIterations', 50, ...
                    'micronsPerDegree', micronsPerDegreeFromOptics, ...
                    'noiseFlag', 'none');
                
            % Set properties
            theOldConeMosaic = setProperties(env, theOldConeMosaic, ...
                varyMPdensity, varyApertureAndOSlength, integrationTime);
        
            % Export
            oldConeMosaicMetaData = theOldConeMosaic.coneData;
            save(fullfile(dataPath, 'coneMosaicHex.mat'), 'theOldConeMosaic', 'oldConeMosaicMetaData');
        else
            load(fullfile(dataPath, 'coneMosaicHex.mat'), 'theOldConeMosaic');
        end
    else
        load(fullfile(dataPath, 'coneMosaicHex.mat'),  'oldConeMosaicMetaData');
        % Generate new mosaic with the old cone mosaic metadata and 300 um/deg
        theNewConeMosaic = cMosaic(...
            'coneData', oldConeMosaicMetaData, ...
            'micronsPerDegree', micronsPerDegreeFromOptics, ...
            'noiseFlag', 'none');
        
        % Set properties
        theNewConeMosaic = setProperties(env, theNewConeMosaic, ...
            varyMPdensity, varyApertureAndOSlength, integrationTime);
        
        % Now generate a new mosaic without the old mosaic metadata
        theNewConeMosaicFull = cMosaic(...
            'eccentricityDegs', [0 0], ...
            'coneDensities', [0.63 0.32 0.05 0.0], ...
            'sizeDegs', mosaicFOVdegs, ...
            'micronsPerDegree', micronsPerDegreeFromOptics, ...
            'noiseFlag', 'none');
        
        % Set properties
        theNewConeMosaicFull = setProperties(env, theNewConeMosaicFull, ...
            varyMPdensity, varyApertureAndOSlength, integrationTime);
    end
    
    
    if (strcmp(env, 'old'))
        % Visualize old mosaic
        hFig = figure(1); clf;
        ax = subplot(1,3,1); cla(ax);
        theOldConeMosaic.visualizeGrid('axesHandle', ax, 'ticksInMicrons', true, ...
                    'visualizedConeAperture',  'geometricArea');
    else
        hFig = figure(1);
        ax = subplot(1,3,2); cla(ax);
        theNewConeMosaic.visualize('figureHandle', hFig, ...
                'axesHandle', ax, 'domain', 'microns', ...
                'domainVisualizationTicks', struct('x', -200:50:200, 'y', -200:50:200), ...
                'crossHairsOnMosaicCenter', true);
            
        ax = subplot(1,3,3); cla(ax);
        theNewConeMosaicFull.visualize('figureHandle', hFig, ...
                'axesHandle', ax, 'domain', 'microns', ...
                'domainVisualizationTicks', struct('x', -200:50:200, 'y', -200:50:200), ...
                'crossHairsOnMosaicCenter', true);
            
        fprintf('cones num (old mosaic): %d (full new mosaic): %d\n', ...
            size(theNewConeMosaic.coneRFpositionsMicrons,1), ...
            size(theNewConeMosaicFull.coneRFpositionsMicrons,1));
    end
    

  
    nTrials = 1; % 30000
    nRepeats = 1;
    noisyResponseInstancesNum = 10;
    
    emPath = zeros(nTrials, 1, 2);
    if (strcmp(env, 'old'))
        tic
        for k = 1:nRepeats
            coneExcitations = theOldConeMosaic.compute(theOI, 'emPath', emPath);
        end
        activationMap = squeeze(coneExcitations(1,:,:));
        fprintf(2,'Old mosaic compute time: %2.2f seconds, max response: %f\n', toc, max(activationMap(:)));
        
        % Also compute 10 noisy response instances
        theOldConeMosaic.noiseFlag = 'random';
        noisyConeExcitationInstances = theOldConeMosaic.compute(theOI, ...
            'emPath', zeros(noisyResponseInstancesNum, 1, 2), ...
            'seed', 12345);
        lConeIndices = find(theOldConeMosaic.pattern == 2);
        mConeIndices = find(theOldConeMosaic.pattern == 3);
        sConeIndices = find(theOldConeMosaic.pattern == 4);
        conesNum = numel(lConeIndices) + numel(mConeIndices) + numel(sConeIndices);
        activationMetaDataOld.noisyConeExcitationInstances = zeros(noisyResponseInstancesNum, conesNum);
        for instanceNo = 1:noisyResponseInstancesNum
            tmp = squeeze(noisyConeExcitationInstances(instanceNo, :,:));
            idx = 1:numel(lConeIndices);
            activationMetaDataOld.noisyConeExcitationInstances(instanceNo,idx) = tmp(lConeIndices);
            activationMetaDataOld.lConeActivations = activationMap(lConeIndices);
            idx = 1:numel(mConeIndices);
            activationMetaDataOld.noisyConeExcitationInstances(instanceNo,numel(lConeIndices)+idx) = tmp(mConeIndices);
            activationMetaDataOld.mConeActivations = activationMap(mConeIndices);
            idx = 1:numel(sConeIndices);
            activationMetaDataOld.noisyConeExcitationInstances(instanceNo,numel(lConeIndices)+numel(mConeIndices)+idx) = tmp(sConeIndices);
            activationMetaDataOld.sConeActivations = activationMap(sConeIndices);
        end

        % Visualize old  mosaic response
        hFig = figure(2);
        ax = subplot(1,3,1);
        cla(ax)
        activationMetaDataOld.meanResponses = theOldConeMosaic.renderActivationMap(ax, activationMap, ...
                'visualizedConeAperture', 'geometricArea', 'mapType', 'modulated disks');
        save(fullfile(dataPath, 'oldMosaicActivationMetaData.mat'), 'activationMetaDataOld');
    else
        tic
        theNewConeMosaic.noiseFlag = 'random';
        for k = 1:nRepeats
            [coneExcitations, noisyConeExcitationInstances] = theNewConeMosaic.compute(theOI, ...
                'seed', 12345, ...
                'nTrials', noisyResponseInstancesNum);
        end

        activationMap = squeeze(coneExcitations(1,1,:));
        fprintf(2,'New mosaic compute time: %2.2f seconds, max response: %f\n', toc, max(activationMap(:)));
        
        
        % Compute with the full new cone mosaic
        tic
        coneExcitationsFull = theNewConeMosaicFull.compute(theOI);
        activationMapFull = squeeze(coneExcitationsFull(1,1,:));
        fprintf(2,'Full mosaic compute time: %2.2f seconds, max response: %f\n', toc, max(activationMapFull(:)));
        
        % Report integrated response over all cones to compare the mosaic
        % coverage
        fprintf(2,'Integrated response: new mosaic (%f), full mosaic (%f)\n', sum(activationMap(:)), sum(activationMapFull(:)));
        
        % Visualize new  mosaic response
        
        figure(2)
        ax = subplot(1,3,2);
        cla(ax)
        theNewConeMosaic.visualize(...
                'figureHandle', hFig, ...
                'axesHandle', ax, ...
                'domain', 'degrees', ...
                'domainVisualizationTicks', struct('x', -1:0.25:1, 'y',  -1:0.25:1), ...
                'activation', activationMap, ...
                'activationRange', [min(activationMap(:)) max(activationMap(:))], ...
                'backgroundColor', [0 0 0], ...
                'fontSize', 18);
            
        ax = subplot(1,3,3);
        cla(ax)
        theNewConeMosaicFull.visualize(...
                'figureHandle', hFig, ...
                'axesHandle', ax, ...
                'domain', 'degrees', ...
                'domainVisualizationTicks', struct('x', -1:0.25:1, 'y', -1:0.25:1), ...
                'activation', activationMapFull, ...
                'activationRange', [min(activationMap(:)) max(activationMap(:))], ...
                'backgroundColor', [0 0 0], ...
                'fontSize', 18);
            
        % Load old activation data
        load(fullfile(dataPath, 'oldMosaicActivationMetaData.mat'), 'activationMetaDataOld');
        
        % Compare activations old vs new mosaic
        activationMetaDataNew.lConeActivations = activationMap(theNewConeMosaic.lConeIndices);
        activationMetaDataNew.mConeActivations = activationMap(theNewConeMosaic.mConeIndices);
        activationMetaDataNew.sConeActivations = activationMap(theNewConeMosaic.sConeIndices);
        
        % Contrast noise response instances
        hFig = figure(4); clf;
        plot(activationMetaDataOld.noisyConeExcitationInstances(:), noisyConeExcitationInstances(:), 'k.');
        xlabel('old mosaic');
        ylabel('new mosaic');
        title('noisy response instances');
        
        % Contrast activations
        hFig = figure(3); clf;
        ax = subplot(2,3,1);
        if (~isempty(activationMetaDataOld.lConeActivations))
            plotActivationCorrespondence(ax, ...
                activationMetaDataOld.lConeActivations(:), ...
                activationMetaDataNew.lConeActivations(:), []);
        end
        
        ax = subplot(2,3,2);
        if (~isempty(activationMetaDataOld.mConeActivations))
            plotActivationCorrespondence(ax, ...
                activationMetaDataOld.mConeActivations(:), ...
                activationMetaDataNew.mConeActivations(:), []);
        end
        
        ax = subplot(2,3,3);
        if (~isempty(activationMetaDataOld.sConeActivations))
            plotActivationCorrespondence(ax, ...
                activationMetaDataOld.sConeActivations(:), ...
                activationMetaDataNew.sConeActivations(:), []);
        end
        
        ax = subplot(2,3,4);
        if (~isempty(activationMetaDataOld.lConeActivations))
            diff = 100*(activationMetaDataNew.lConeActivations - activationMetaDataOld.lConeActivations)./ activationMetaDataNew.lConeActivations;
            plotActivationCorrespondence(ax, ...
                activationMetaDataOld.lConeActivations(:), ...
                diff(:), ...
                3*[-1 1]);
        end
        
        ax = subplot(2,3,5);
        if (~isempty(activationMetaDataOld.mConeActivations))
            diff = 100*(activationMetaDataNew.mConeActivations - activationMetaDataOld.mConeActivations)./ activationMetaDataNew.mConeActivations;
            plotActivationCorrespondence(ax, ...
                activationMetaDataOld.mConeActivations(:), ...
                diff(:), ...
                3*[-1 1]);
        end
        
        ax = subplot(2,3,6);
        if (~isempty(activationMetaDataOld.sConeActivations))
            diff = 100*(activationMetaDataNew.sConeActivations - activationMetaDataOld.sConeActivations)./ activationMetaDataNew.sConeActivations;
            plotActivationCorrespondence(ax, ...
                activationMetaDataOld.sConeActivations(:), ...
                diff(:), ...
               3*[-1 1]);
        end
        
       NicePlot.exportFigToPDF('residuals.pdf', hFig, 300);
    end
    
end

function  theConeMosaic = setProperties(env, theConeMosaic, ...
        varyMPdensity, varyApertureAndOSlength, integrationTime)
    
    % Integration time
    theConeMosaic.integrationTime = integrationTime;
    
    % MP density
    if (strcmp(env, 'old'))
        theConeMosaic.apertureBlur = true;
        theConeMosaic.eccBasedMacularPigment = varyMPdensity;
    else
        theConeMosaic.eccVaryingMacularPigmentDensity = varyMPdensity;
        theConeMosaic.eccVaryingMacularPigmentDensityDynamic = false;
    end
    
    
    % aperture and OS length
    if (strcmp(env, 'old'))
        theConeMosaic.eccBasedConeQuantalEfficiency = varyApertureAndOSlength;
    else
        theConeMosaic.eccVaryingConeAperture = varyApertureAndOSlength;
        theConeMosaic.eccVaryingOuterSegmentLength = varyApertureAndOSlength;
        theConeMosaic.eccVaryingConeBlur = false;
    end

    
    
end

function plotActivationCorrespondence(ax, oldActivations, newActivations, yRange)
        
    emptyYRange = false;
    if (isempty(yRange)) 
        emptyYRange = true;
        xxL(1) = min([min(oldActivations) min(newActivations)]);
        xxL(2) = max([max(oldActivations) max(newActivations)]);
        if (xxL(1) == xxL(2))
            xxL(1) = xxL(2)-1;
            xxL(2) = xxL(2)+1;
        end
        yRange = xxL;
    else
        xxL(1) = min(oldActivations);
        xxL(2) = max(oldActivations);
        if (xxL(1) == xxL(2))
            xxL(1) = xxL(2)-1;
            xxL(2) = xxL(2)+1;
        end
    end
    
    if (emptyYRange)  
        plot(ax,[0 10000], [0 10000], 'r-', 'LineWidth', 1.0);
    else
        plot(ax,[0 10000], [0 0], 'r-', 'LineWidth', 1.0);
    end

    hold(ax, 'on');
    plot(ax,oldActivations, newActivations, 'k.');
    xlabel(ax,'absorptions (@coneMosaicHex)')
    
    if (emptyYRange)   
        ylabel(ax,'absorptions (@cMosaic)');
    else
        ylabel(ax, 'error');
    end
    set(ax, 'XLim', xxL, 'YLim', yRange);
    grid(ax, 'on')
    axis(ax, 'square');
    set(ax, 'FontSize', 12)
end


    
function [theOI, micronsPerDegreeFromOptics] = computeTheOI(dataPath)
    
    % Load scene
    d = load(fullfile(dataPath, 'scene10cpd_15_2_5degs_optimalPixelSize.mat'), 'theScene');
    theScene = d.theScene;
    
    fprintf('Computing OI\n');
    theOI = oiCreate('wvf human');
    theOI = oiCompute(theOI, theScene);
    
    oiSize = oiGet(theOI, 'size');
    % Compute optics micronsPerDegree
    horizontalFOVdegs = oiGet(theOI, 'wangular');
    oiResMicrons = oiGet(theOI, 'height spatial resolution')*1e6;
    oiResDegs = horizontalFOVdegs/oiSize(2);
    micronsPerDegreeFromOptics = oiResMicrons/oiResDegs;
            
    fprintf('Done computing OI\n');
    
end

