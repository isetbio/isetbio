function coneMosaic(app)

    % Generate cone mosaic with new params
    app.components.coneMosaic = cMosaic(...
        'whichEye', app.coneMosaicParams.whichEye, ...
        'sizeDegs', app.coneMosaicParams.sizeDegs, ...
        'eccentricityDegs', app.coneMosaicParams.eccentricityDegs);
    
    conePos = bsxfun(@minus, app.components.coneMosaic.coneRFpositionsDegs, app.components.coneMosaic.eccentricityDegs);
    coneEccArcMin = 60*sqrt(sum(conePos.^2,2));
    [~, idx] = sort(coneEccArcMin);
    centralCone = idx(1);
    conePos = bsxfun(@minus, app.components.coneMosaic.coneRFpositionsDegs, app.components.coneMosaic.coneRFpositionsDegs(centralCone,:));
    coneEccArcMin = 60*sqrt(sum(conePos.^2,2));
    [~, idx] = sort(coneEccArcMin);
    
    % Extract the centralConeOutlinesArcMin 
    for k = 1:app.centralConeOutlinesNum
        coneIndex = idx(k);
        radius = 0.5*app.components.coneMosaic.coneRFspacingsDegs(coneIndex)*60*app.components.coneMosaic.coneApertureToDiameterRatio;
        app.centralConeOutlinesArcMin(k,1,:) = (conePos(coneIndex,1)-conePos(centralCone,1))*60 + radius * cosd(0:5:360);
        app.centralConeOutlinesArcMin(k,2,:) = (conePos(coneIndex,2)-conePos(centralCone,2))*60 + radius * sind(0:5:360);
    end
    
    
    
%     , ...
%         'coneDensities', app.mosaicConeDensities, ...
%         'tritanopicRadiusDegs', app.mosaicTritanopicRadiusDegs, ...
%         'wave', app.wavelengthMinValue : app.wavelengthStepSize: app.wavelengthMaxValue, ...
%         'integrationTime', app.mosaicIntegrationTimeMseconds/1000, ...
%         'eccVaryingConeAperture', app.mosaicEccVaryingConeAperture, ...
%         'eccVaryingOuterSegmentLength', app.mosaicEccVaryingOuterSegmentLength, ...
%         'eccVaryingConeBlur', app.mosaicEccVaryingConeBlur, ...
%         'eccVaryingMacularPigmentDensity', app.mosaicEccVaryingMacularPigmentDensity, ...
%         'eccVaryingMacularPigmentDensityDynamic', false ...


    %   Visualize the cone mosaic
    CSFGeneratorApp.render.coneMosaicView(app);
    

    %   Set the microns-per-deg field
    app.visualFieldMagnificationFactorEditField.Value = app.components.coneMosaic.micronsPerDegree;   
                

    
end

