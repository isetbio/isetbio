function [theStimulusFxCoord, theStimulusFyCoord, stimulusInInsideAnalyzedSFregion, theStimulusFrequency, spatialFrequencySupport] = ...
    HartleySFmap(HartleySpatialModulationPatterns, spatialSupportDegs, lIndices, mIndices, theStimulusIndex, ...
    visualizeStimulusSet, maxSFtoBeAnalyzed)

    nStimuli = numel(lIndices);
    omega = max(lIndices);
    assert(nStimuli == (2*omega+1)^2, 'error in decoding Hartley stimuli');

    stimSizeDegs = max(spatialSupportDegs)-min(spatialSupportDegs);
    sf = zeros(1, numel(lIndices));
    spatialFrequencySupportXgrid = zeros(2*omega+1, 2*omega+1);
    spatialFrequencySupportYgrid = zeros(2*omega+1, 2*omega+1);

    ixList = zeros(1, numel(lIndices));
    iyList = ixList;

    for sIndex = 1:numel(lIndices)
        fx = lIndices(sIndex)/stimSizeDegs;
        fy = mIndices(sIndex)/stimSizeDegs;
        
        ix = lIndices(sIndex)+omega+1;
        iy = mIndices(sIndex)+omega+1;
       
        ixList(sIndex) = ix;
        iyList(sIndex) = iy;
        spatialFrequencySupportXgrid(iy,ix) = fx;
        spatialFrequencySupportYgrid(iy,ix) = fy;
        sf(sIndex) = sqrt(fx^2+fy^2);

        if (theStimulusIndex == sIndex)
            if (abs(sf(sIndex)) <= maxSFtoBeAnalyzed)
                stimulusInInsideAnalyzedSFregion = true;
            else
                stimulusInInsideAnalyzedSFregion = false;
            end
            theStimulusFxCoord = ix;
            theStimulusFyCoord = iy;
            theStimulusFrequency = sf(sIndex);
        end
    end

    
    spatialFrequencySupport = spatialFrequencySupportXgrid(1,:);

    stimuliIndicesKept = find(sf <= maxSFtoBeAnalyzed);
    if (visualizeStimulusSet) && (numel(stimuliIndicesKept)<25*25)

        ixListKept = ixList(stimuliIndicesKept);
        iyListKept = iyList(stimuliIndicesKept);
        
        rowsNum = max(iyListKept)-min(iyListKept)+1;
        colsNum = max(ixListKept)-min(ixListKept)+1;
        hFig = figure(1001); clf;
        set(hFig, 'Position', [10 10 1280 1350], 'Color', [1 1 1]);
        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', rowsNum, ...
           'colsNum', colsNum, ...
           'heightMargin',  0.001, ...
           'widthMargin',    0.001, ...
           'leftMargin',     0.0, ...
           'rightMargin',    0.0, ...
           'bottomMargin',   0.0, ...
           'topMargin',      0.0);

        for ii = 1:numel(stimuliIndicesKept)
            sIndex = stimuliIndicesKept(ii);

            subplotCol = ixListKept(ii)-min(ixListKept)+1;
            subplotRow = iyListKept(ii)-min(iyListKept)+1;
            

            ax = subplot('Position', subplotPosVectors(subplotRow, subplotCol).v);
            imagesc(ax,spatialSupportDegs, spatialSupportDegs, squeeze(HartleySpatialModulationPatterns(sIndex, :,:)));
            set(ax, 'XTick', [], 'YTick', []);
            set(ax, 'CLim', [-1 1]);
            axis(ax, 'image')
            colormap(ax, brewermap(1024, '*greys'));
            title(ax, sprintf('%2.1fc/deg', sf(sIndex)))
            drawnow;
        end

    end

    visualizeFullStimulusSet = false;
    if (visualizeFullStimulusSet)

        hFig = figure(1001); clf;
        set(hFig, 'Position', [10 10 1280 1350], 'Color', [1 1 1]);
        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 2*omega+1, ...
           'colsNum', 2*omega+1, ...
           'heightMargin',  0.001, ...
           'widthMargin',    0.001, ...
           'leftMargin',     0.0, ...
           'rightMargin',    0.0, ...
           'bottomMargin',   0.0, ...
           'topMargin',      0.0);
    

        for sIndex = 1:numel(lIndices)
            subplotCol = lIndices(sIndex)+omega+1;
            subplotRow = mIndices(sIndex)+omega+1;

            ax = subplot('Position', subplotPosVectors(subplotRow, subplotCol).v);
            imagesc(ax,spatialSupportDegs, spatialSupportDegs, squeeze(HartleySpatialModulationPatterns(sIndex, :,:)));
            
            if (abs(sf(sIndex)) <= maxSFtoBeAnalyzed)
                hold(ax, 'on')
                xx = [spatialSupportDegs(1) spatialSupportDegs(1) spatialSupportDegs(end) spatialSupportDegs(end) spatialSupportDegs(1)];
                yy = [spatialSupportDegs(1) spatialSupportDegs(end) spatialSupportDegs(end) spatialSupportDegs(1) spatialSupportDegs(1)];
                plot(ax, xx,yy, 'r-');
                hold(ax, 'off')
                if (abs(sf(sIndex))>10)
                    title(ax, sprintf('%2.0fc/deg', sf(sIndex)))
                else
                    title(ax, sprintf('%2.1fc/deg', sf(sIndex)))
                end

            end
            set(ax, 'XTick', [], 'YTick', []);
            set(ax, 'CLim', [-1 1]);
            axis(ax, 'image')
            colormap(ax, brewermap(1024, '*greys'));
            drawnow;
        end
    end


end

