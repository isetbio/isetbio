function visualizeOpticsAtEccentricities(obj, eccDegs, opticsParams)

    assert(size(eccDegs,1)<=9, 'The number of visualized positions must be <= 9');

    % Visualized PSF range and wavelength
    psfRangeDegs = 5/60;
    targetWavelength = 550;
    [~,idx] = min(abs(obj.inputConeMosaic.wave-targetWavelength));

    % Plot format
    ff = MSreadyPlot.figureFormat('3x3');
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 ff.figureSize(1) ff.figureSize(2)], 'Color', [1 1 1]);
    theAxes = MSreadyPlot.generateAxes(hFig,ff);

    % Plot rows and cols
    [X,Y] = meshgrid(1:3, 1:3); cols = X(:); rows = Y(:);

    % Render PSF at each position
    for iPos = 1:size(eccDegs,1)
        % Generate optics at the current position
        opticsParams.positionDegs = eccDegs(iPos,:);
        dataOut = obj.generateOptics(opticsParams);

        % Retrieve the PSF data
        thePSFData = dataOut.thePSFData;
        thePSFData.psfSupportXdegs = thePSFData.supportX/60;
        thePSFData.psfSupportYdegs = thePSFData.supportY/60;
        thePSFData.data = squeeze(thePSFData.data(:,:,idx));

        % Render the PSF
        plotTitle = sprintf('XYecc (degs): (%2.1f, %2.1f)', ...
            opticsParams.positionDegs(1), opticsParams.positionDegs(2));
        noXLabel = (rows(iPos) < 3);
        noYLabel = (cols(iPos)> 1);
        MSreadyPlot.render2DPSF(theAxes{rows(iPos),cols(iPos)}, ...
            thePSFData.psfSupportXdegs, thePSFData.psfSupportYdegs, ...
            thePSFData.data, psfRangeDegs, plotTitle, ff, ...
            'noXLabel', noXLabel, ...
            'noYLabel', noYLabel);
        drawnow;
    end

end



