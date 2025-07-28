function hFig = visualizeOpticsAtEccentricities(obj, eccDegs, opticsParams, tickSeparationArcMin)

    assert(size(eccDegs,1)<=9, 'The number of visualized positions must be <= 9');

    % Visualized PSF range and wavelength
    psfRangeDegs = 0.5*(tickSeparationArcMin*4)/60;
    targetWavelength = 550;
    [~,idx] = min(abs(obj.inputConeMosaic.wave-targetWavelength));


    % Plot format
    ff = MSreadyPlot.figureFormat('3x3');
    hFig = figure(1); clf;
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);
    
    % Plot rows and cols
    [X,Y] = meshgrid(1:3, 1:3); 
    cols = X(:); rows = Y(:);
    rows = 4 - rows;

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

        % Flip left-right
        thePSFData.data = fliplr(thePSFData.data);

        % Render the PSF
        plotTitle = sprintf('XYecc (degs): (%2.1f, %2.1f)', ...
            opticsParams.positionDegs(1), opticsParams.positionDegs(2));
        noXLabel = (rows(iPos) < 3);
        noYLabel = (cols(iPos)> 1);

        apertureDataStruct = localConeApertureData(obj, opticsParams.positionDegs, thePSFData.psfSupportXdegs);

        MSreadyPlot.render2DPSF(theAxes{rows(iPos),cols(iPos)}, ...
            thePSFData.psfSupportXdegs, thePSFData.psfSupportYdegs, ...
            thePSFData.data, psfRangeDegs, plotTitle, ff, ...
            'withConeApertureData', apertureDataStruct, ...
            'tickSeparationArcMin', tickSeparationArcMin, ...
            'noXLabel', noXLabel, ...
            'noYLabel', noYLabel);
        drawnow;
    end

end


function dOut = localConeApertureData(obj, opticsPositionDegs, psfSupportDegs)
    d = sqrt(sum((bsxfun(@minus, obj.inputConeMosaic.coneRFpositionsDegs, opticsPositionDegs)).^2,2));
    idx = find (d<max(psfSupportDegs));
    coneAperturePositionsDegs = bsxfun(@minus, obj.inputConeMosaic.coneRFpositionsDegs(idx,:),opticsPositionDegs);

    if (isfield(obj.inputConeMosaic.coneApertureModifiers, 'shape') && (strcmp(obj.inputConeMosaic.coneApertureModifiers.shape, 'Gaussian')))
        coneAperturesDegs = obj.inputConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor * ...
                            obj.inputConeMosaic.coneApertureDiametersDegs(idx);
    else  
        coneAperturesDegs = obj.inputConeMosaic.coneApertureDiametersDegs(idx);
    end

    dOut = struct(...
       'positionDegs', coneAperturePositionsDegs, ...
       'RcDegs', coneAperturesDegs);

end

