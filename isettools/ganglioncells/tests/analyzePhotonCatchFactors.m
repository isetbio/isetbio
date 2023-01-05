function hFigs = analyzePhotonCatchFactors()
    % Load ISETBio data across eccentricities
    load('ISETBioData.mat', 'dataDict', 'eccDegsGrid');

    % Retrieve the metadata
    dMetaDataStruct = dataDict('metaData');

    % Retrieve the foveal data
    dataLabel = sprintf('eccXY = %2.2f,%2.2f', 0, 0);
    dStruct = dataDict(dataLabel);

    % Foveal values (reference)
    fovealInnerSegmentDiameterDegs = dStruct.innerSegmentDiameterDegs;
    fovealOuterSegmentLengthMicrons = dStruct.outerSegmentLengthMicrons;
    fovealAbsorptanceSpectra = dStruct.absorptanceSpectra;
    fovealMacularPigmentTransmittance = dStruct.macularPigmentTransmittance;

    % Retrieve data at all eccentricities
    for iEcc = 1:size(eccDegsGrid,1)

        % Retrieve data for this eccentricity
        dataLabel = sprintf('eccXY = %2.2f,%2.2f', eccDegsGrid(iEcc,1), eccDegsGrid(iEcc,2));
        dStruct = dataDict(dataLabel);

        % How inner segment diameter affects photon catch rate (relative to the foveal cone )
        innerSegmentDiameterCatchFactor(iEcc) = (dStruct.innerSegmentDiameterDegs/fovealInnerSegmentDiameterDegs).^2;

        % How outer segment length affects photon catch rate (relative to the foveal cone)
        outerSegmentLengthCatchFactors(iEcc) = dStruct.outerSegmentLengthMicrons  / fovealOuterSegmentLengthMicrons;

        % How absorptance, integrated over wavelength, affects photon catch rate (relative to the
        % foveal cone), separately for each cone type
        absorptanceCatchFactors(iEcc,:) = sum(dStruct.absorptanceSpectra,2) ./ sum(fovealAbsorptanceSpectra,2);

        % How macular pigment, integrated over wavelength) affects photon catch rate (relative to the
        % foveal cone)
        macularPigmentTransmittanceFactors(iEcc) = sum(dStruct.macularPigmentTransmittance,2) / sum(fovealMacularPigmentTransmittance,2);

        % How absorptance * macular pigment transmittance (integrated over wavelength) affects photon catch rate, 
        % separately for each cone type
        for coneTypeIndex = cMosaic.LCONE_ID:cMosaic.SCONE_ID
            fovealPhotonCatch = fovealAbsorptanceSpectra(coneTypeIndex,:) .* fovealMacularPigmentTransmittance;
            photonCatchAtTargetEccentricity = dStruct.absorptanceSpectra(coneTypeIndex,:) .* dStruct.macularPigmentTransmittance;
            % Integrate over wavelength
            photonCatchFactors(iEcc,coneTypeIndex) = sum(photonCatchAtTargetEccentricity) / sum(fovealPhotonCatch);
        end

        % Absortance * macular pigment transmittance * inner segment diameter (total photon catch factor), 
        % separately for the 3 cone types
        singleConePhotonCatchFactors(iEcc,:) = photonCatchFactors(iEcc,:) * innerSegmentDiameterCatchFactor(iEcc);
        
        % Total photon catch factor in mRGC RF center (only L/M cones in RF center)
        midgetRGCRFcenterConeTypes = [dMetaDataStruct.LCONE_ID dMetaDataStruct.MCONE_ID];
        meanLMconePhotonCatchFactor = mean(singleConePhotonCatchFactors(iEcc,midgetRGCRFcenterConeTypes),2);
        mRGCRFcenterCatchFactors(iEcc) = meanLMconePhotonCatchFactor * dStruct.meanConesNumInRFcenter;
    end % iEcc


    % Plot the variation of different factors along the horizontal meridian
    hFig1 = figure(1); clf;
    set(hFig1, 'Position', [10 10 1500 740], 'Color', [1 1 1]);
    
    % Find all eccentricities with eccY = 0
    idx = find(eccDegsGrid(:,2) == 0);
    eccXdegs = eccDegsGrid(idx,1);

    ax = subplot(1,2,1);
    hold (ax, 'on');
    plot(ax, eccXdegs, log10(absorptanceCatchFactors(idx,dMetaDataStruct.LCONE_ID)), ...
        'ro-', 'LineWidth', 1.5, 'MarkerSize', 14, 'MarkerFaceColor', [1 0.5 0.5]); 
    plot(ax, eccXdegs, log10(absorptanceCatchFactors(idx,dMetaDataStruct.MCONE_ID)), ...
        'go-', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', [0.5 0.8 0.5]); 
    plot(ax, eccXdegs, log10(absorptanceCatchFactors(idx,dMetaDataStruct.SCONE_ID)), ...
        'bo-', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', [0.5 0.5 1]); 
    plot(ax, eccXdegs, log10(innerSegmentDiameterCatchFactor(idx)), ...
        'ko-', 'LineWidth', 1.5, 'MarkerSize', 12, 'MarkerFaceColor', [0.7 0.7 0.7]); 
    plot(ax, eccXdegs, log10(macularPigmentTransmittanceFactors(idx)), ...
        'mo-', 'LineWidth', 1.5, 'MarkerSize', 12, 'MarkerFaceColor', [0.7 0.5 0.7]); 
    
    set(ax, 'YLim', [-0.3 1.5], 'YTick', -1:0.1:2, 'FontSize', 16, 'TickDir', 'both');
    grid(ax, 'on'); box(ax, 'off');
    legend(ax, {'L-cone absorptance', ...
                'M-cone absorptance', ...
                'S-cone absorptance', ...
                'inner segment diameter', ...
                'macular pigment transmittance' ...
                }, ...
                'Location', 'North', 'Box', 'off');
    ylabel(ax, 'log10(relative catch rate)');
    xlabel(ax, 'eccentricity, x (degs)')

    ax = subplot(1,2,2);
    plot(ax, eccXdegs, log10(singleConePhotonCatchFactors(idx,dMetaDataStruct.LCONE_ID)), ...
        'ro-', 'LineWidth', 1.5, 'MarkerSize', 14, 'MarkerFaceColor', [1 0.5 0.5]); 
    hold (ax, 'on');
    plot(ax, eccXdegs, log10(singleConePhotonCatchFactors(idx,dMetaDataStruct.MCONE_ID)), ...
        'go-', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', [0.5 0.8 0.5]); 
    plot(ax, eccXdegs, log10(singleConePhotonCatchFactors(idx,dMetaDataStruct.SCONE_ID)), ...
        'bo-', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', [0.5 0.5 1]); 
    plot(ax, eccXdegs, log10(mRGCRFcenterCatchFactors(idx)), ...
        'ko-', 'LineWidth', 1.5, 'MarkerSize', 12, 'MarkerFaceColor', [0.7 0.7 0.7]); 
    set(ax, 'YLim', [-0.3 1.5], 'YTick', -1:0.1:2, 'FontSize', 16,  'TickDir', 'both');
    grid(ax, 'on'); box(ax, 'off');
    legend(ax, {'L-cone catch (absorptance + inner segment diam + mac. pigm.)', ...
                'M-cone catch (absorptance + inner segment diam + mac. pigm.)', ...
                'S-cone catch (absorptance + inner segment diam + mac. pigm.)', ...
                'mRGC RFcenter catch'}, ...
                'Location', 'North', 'Box', 'off');
    xlabel(ax, 'eccentricity, x (degs)')


    % Interpolate data
    x = eccDegsGrid(:,1);
    y = eccDegsGrid(:,2);
    intepolationMethod = 'natural';  % 'linear'
    xInterp = min(x):0.2:max(x);
    yInterp = min(y):0.2:max(y);
    [X,Y] = meshgrid(xInterp, yInterp);
   
    singleLConeCatchFactorMap = griddata(x,y,singleConePhotonCatchFactors(:,dMetaDataStruct.LCONE_ID)',X,Y, intepolationMethod);
    singleMConeCatchFactorMap = griddata(x,y,singleConePhotonCatchFactors(:,dMetaDataStruct.MCONE_ID)',X,Y, intepolationMethod);
    singleSConeCatchFactorMap = griddata(x,y,singleConePhotonCatchFactors(:,dMetaDataStruct.SCONE_ID)',X,Y, intepolationMethod);
    mRGCRFcenterCatchFactorMap = griddata(x,y,mRGCRFcenterCatchFactors',X,Y, intepolationMethod);


    singleConeRelativeCatchRateRange  = [0 ceil(max(singleConePhotonCatchFactors(:)))];
    rfCenterRelativeCatchRateRange = [0 ceil(max(mRGCRFcenterCatchFactorMap(:)))];

    zLevelsSingleCone = singleConeRelativeCatchRateRange(1):0.25:singleConeRelativeCatchRateRange(2);
    zLevelsMidgetRFcenter = rfCenterRelativeCatchRateRange(1):0.5:rfCenterRelativeCatchRateRange(2);

    rfCenterRelativeCatchRateRange = [1 25];

  
    hFig2 = figure(2); clf;
    set(hFig2, 'Position', [10 10 1450 1000], 'Color', [1 1 1]);
    
    ax = subplot(2,2,1);
    plotCatchFactorMap(ax, X,Y,singleLConeCatchFactorMap, zLevelsMidgetRFcenter, ...
        rfCenterRelativeCatchRateRange, dMetaDataStruct.cMap, false, true, 'L-cone catch');


    ax = subplot(2,2,2);
    plotCatchFactorMap(ax, X,Y,singleMConeCatchFactorMap, zLevelsMidgetRFcenter, ...
        rfCenterRelativeCatchRateRange, dMetaDataStruct.cMap, false, false, 'M-cone catch');

    ax = subplot(2,2,3);
    plotCatchFactorMap(ax, X,Y,singleSConeCatchFactorMap, zLevelsMidgetRFcenter, ...
        rfCenterRelativeCatchRateRange, dMetaDataStruct.cMap, true, true, 'S-cone catch');

    ax = subplot(2,2,4);
    plotCatchFactorMap(ax, X,Y, mRGCRFcenterCatchFactorMap, zLevelsMidgetRFcenter, ...
        rfCenterRelativeCatchRateRange, dMetaDataStruct.cMap, true, false, 'mRGC RFcenter catch');

    hFigs = {hFig1, hFig2};
end

function plotCatchFactorMap(ax, X,Y, catchFactorMap, zLevels, catchRateRange, cMap, labelXmeridian, labelYmeridian, plotTitle)

    contourf(ax,X,Y,catchFactorMap, zLevels);
    c = colorbar;
    c.Label.String = 'relative catch rate ';
    axis(ax, 'equal');
    grid(ax, 'on');
    box(ax, 'off');
    set(ax, 'XLim', [-8 8], 'YLim', [-6 6], 'XTick', -10:1:10, 'YTick', -10:1:10);
    set(ax, 'XColor', [0.3 0.3 0.3], 'YColor', [0.3 0.3 0.3], 'LineWidth', 1.0, 'CLim', catchRateRange);
    set(ax, 'TickDir', 'both');
    set(ax, 'FontSize', 16);
    if (labelXmeridian)
        xlabel(ax, '\leftarrow  {\it temporal retina}        (degs)        {\it  nasal retina} \rightarrow       ');
    end
    if (labelYmeridian)
        ylabel(ax,'\leftarrow  {\it inferior retina}  (degs) {\it superior retina} \rightarrow ');
    end
    
    colormap(cMap);
    xtickangle(ax, 0);
    title(ax, plotTitle)
end
