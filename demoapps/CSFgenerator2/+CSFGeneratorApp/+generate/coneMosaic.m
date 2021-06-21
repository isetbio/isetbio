function theConeMosaic = coneMosaic(app, dialog)

    if (isa(app, 'ISETBioCSFGenerator'))
        appCall = true;
    else
        appCall = false;
    end
    
    if (appCall)
        deleteProgressBar = isempty(dialog);
        if (deleteProgressBar)
            % Open progressbar
            dialogBox = uiprogressdlg(app.mainView,'Title','Please Wait',...
                        'Message','Generating cone mosaic ...');
            dialogBox.Value = 0.2; 
        end
    end
    
    
    % Generate cone mosaic with new params
    coneDensities = [app.coneMosaicParams.lConeRatio  app.coneMosaicParams.mConeRatio  app.coneMosaicParams.sConeRatio];
    wavelengths = app.stimParams.wavelengthSupportMin:app.stimParams.wavelengthSupportStepSize:app.stimParams.wavelengthSupportMax;
    app.components.coneMosaic = cMosaic(...
        'wave', wavelengths, ...
        'whichEye', app.coneMosaicParams.whichEye, ...
        'sizeDegs', app.coneMosaicParams.sizeDegs, ...
        'eccentricityDegs', app.coneMosaicParams.eccentricityDegs, ...
        'coneDensities', coneDensities/sum(coneDensities), ...
        'tritanopicRadiusDegs', app.coneMosaicParams.tritanopicRadiusDegs, ...
        'integrationTime', app.coneMosaicParams.integrationTimeSeconds, ...
        'eccVaryingConeAperture', app.coneMosaicParams.eccVaryingConeAperture, ...
        'eccVaryingOuterSegmentLength', app.coneMosaicParams.eccVaryingOuterSegmentLength, ...
        'eccVaryingConeBlur', app.coneMosaicParams.eccVaryingConeApertureBlur, ...
        'eccVaryingMacularPigmentDensity', app.coneMosaicParams.eccVaryingMacularPigmentDensity, ...
        'eccVaryingMacularPigmentDensityDynamic', app.coneMosaicParams.eccVaryingMacularPigmentDynamic);
    
    if (~appCall)
        theConeMosaic = app.components.coneMosaic;
        return;
    end
    
    % Reset view limits
    app.coneMosaicViewXLimsDegs = app.coneMosaicParams.eccentricityDegs(1) + 0.55*app.coneMosaicParams.sizeDegs(1)*[-1 1];
    app.coneMosaicViewYLimsDegs = app.coneMosaicParams.eccentricityDegs(2) + 0.55*app.coneMosaicParams.sizeDegs(2)*[-1 1];
    
    if (~isempty(app.components.coneMosaic.micronsPerDegreeApproximation))
       app.coneMosaicViewXLimsMicrons = app.coneMosaicViewXLimsDegs * app.components.coneMosaic.micronsPerDegreeApproximation; 
       app.coneMosaicViewYLimsMicrons = app.coneMosaicViewYLimsDegs * app.components.coneMosaic.micronsPerDegreeApproximation;
    else
        for k = 1:2
            app.coneMosaicViewXLimsMicrons(k) = 1e3 * RGCmodels.Watson.convert.rhoDegsToMMs(app.coneMosaicViewXLimsDegs(k));
            app.coneMosaicViewYLimsMicrons(k) = 1e3 * RGCmodels.Watson.convert.rhoDegsToMMs(app.coneMosaicViewYLimsDegs(k));
        end
    end
            
    % Generate outlines for central cones (used when plotting the PSF)
    if (~isempty(app.components.coneMosaic.coneRFpositionsDegs))
        conePos = bsxfun(@minus, app.components.coneMosaic.coneRFpositionsDegs, app.components.coneMosaic.eccentricityDegs);
        coneEccArcMin = 60*sqrt(sum(conePos.^2,2));
        [~, idx] = sort(coneEccArcMin);
        centralCone = idx(1);
        conePos = bsxfun(@minus, app.components.coneMosaic.coneRFpositionsDegs, app.components.coneMosaic.coneRFpositionsDegs(centralCone,:));
        coneEccArcMin = 60*sqrt(sum(conePos.^2,2));
        [~, idx] = sort(coneEccArcMin);

        % Extract the centralConeOutlinesArcMin 
        for k = 1:app.centralConeOutlinesNum
            if (k <= numel(idx))
                coneIndex = idx(k);
                radius = 0.5*app.components.coneMosaic.coneRFspacingsDegs(coneIndex)*60*app.components.coneMosaic.coneApertureToDiameterRatio;
                app.centralConeOutlinesArcMin(k,1,:) = (conePos(coneIndex,1)-conePos(centralCone,1))*60 + radius * cosd(0:5:360);
                app.centralConeOutlinesArcMin(k,2,:) = (conePos(coneIndex,2)-conePos(centralCone,2))*60 + radius * sind(0:5:360);
            else
                app.centralConeOutlinesArcMin(k,1,:) = nan;
                app.centralConeOutlinesArcMin(k,2,:) = nan;
            end
        end
    else
        for k = 1:app.centralConeOutlinesNum
            app.centralConeOutlinesArcMin(k,1,:) = nan;
            app.centralConeOutlinesArcMin(k,2,:) = nan;
        end
    end
        
    % Set the microns-per-deg field
    app.roiMagnificationFactorEditField.Value = app.components.coneMosaic.micronsPerDegree;   

    % Update status for 'cone mosaic'
    conesNum = size(app.components.coneMosaic.coneRFpositionsDegs,1);
    app.statusMessages('cone mosaic')  = struct(...
        'text', sprintf('Generated cone mosaic with %d cones.', conesNum), ...
        'fontColor', app.colors('good message foreground'), ...
        'backgroundColor', app.colors('good message background'), ...
        'fontWeight', 'normal');
    
    % Render the status on the status field of tab B
    CSFGeneratorApp.render.statusField(app,'B', 'cone mosaic'); 
    
    if (deleteProgressBar)
        % Close progressbar
        close(dialogBox);
    end
end

