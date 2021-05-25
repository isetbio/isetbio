function visualizeFullAbsorptionsDensity(obj, figNo)

        hFig = figure(figNo); clf;
        set(hFig, 'Position', [1000 10 1000 1000]);
        
        colormap(gray(1024));
        
        deltaAngle = 15;
        iTheta = (0:deltaAngle:360) / 180 * pi;
        coneApertureShape.x = cos(iTheta);
        coneApertureShape.y = sin(iTheta);
    
        visualizedConeAperture = 'light collecting area';
        if (strcmp(visualizedConeAperture, 'geometricArea'))
            visualizedApertureMultiplier = 1.0;
        else
            visualizedApertureMultiplier = obj.coneApertureToDiameterRatio;
        end
    
        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
         'rowsNum', 2, ...
         'colsNum', 2, ...
         'heightMargin',  0.03, ...
         'widthMargin',    0.03, ...
         'leftMargin',     0.03, ...
         'rightMargin',    0.01, ...
         'bottomMargin',   0.03, ...
         'topMargin',      0.03);
     
        for coneTypeIndex = 1:size(obj.absorptionsDensityFullMap,3)
            
            row = floor((coneTypeIndex-1)/2)+1;
            col = mod((coneTypeIndex-1),2)+1;
            ax = subplot('Position', subplotPosVectors(row,col).v);
            
            imagesc(ax,obj.absorptionsDensitySpatialSupportMicrons{2}, ...
                    obj.absorptionsDensitySpatialSupportMicrons{1}, ...
                    squeeze(obj.absorptionsDensityFullMap(:,:,coneTypeIndex)));
            hold(ax, 'on');
            switch(coneTypeIndex)
                case 1
                    % Plot L-cones
                    renderPatchArray(ax, coneApertureShape, visualizedApertureMultiplier*obj.coneRFspacingsMicrons(obj.lConeIndices)*0.5, ...
                        obj.coneRFpositionsMicrons(obj.lConeIndices,:), [1 0 0], 1.0);
                    
                case 2 
                    % Plot M-cones
                    renderPatchArray(ax, coneApertureShape, visualizedApertureMultiplier*obj.coneRFspacingsMicrons(obj.mConeIndices)*0.5, ...
                        obj.coneRFpositionsMicrons(obj.mConeIndices,:), [0 1 0], 1.0);
                case 3
                    % Plot S-cones
                    renderPatchArray(ax, coneApertureShape, visualizedApertureMultiplier*obj.coneRFspacingsMicrons(obj.sConeIndices)*0.5, ...
                    obj.coneRFpositionsMicrons(obj.sConeIndices,:),  [0 0.5 1], 1.0);
            end
        
            hold(ax, 'off');
            axis(ax,'xy');
            axis(ax,'equal');
            minXY = min(obj.coneRFpositionsMicrons,[],1);
            maxXY = max(obj.coneRFpositionsMicrons,[],1);
            set(ax, 'XLim', [minXY(1) maxXY(1)], 'YLim', [minXY(2) maxXY(2)]);
        end

end


function renderPatchArray(axesHandle, apertureShape, apertureRadii, rfCoords, ...
    edgeColor, lineWidth)

    conesNum = numel(apertureRadii);
    if (conesNum == 0)
        return;
    end
    
    verticesPerCone = numel(apertureShape.x);
    
    verticesList = zeros(verticesPerCone * conesNum, 2);
    facesList = [];
    
    for coneIndex = 1:conesNum
        idx = (coneIndex - 1) * verticesPerCone + (1:verticesPerCone);
        verticesList(idx, 1) = apertureShape.x*apertureRadii(coneIndex) + rfCoords(coneIndex,1);
        verticesList(idx, 2) = apertureShape.y*apertureRadii(coneIndex) + rfCoords(coneIndex,2);

        facesList = cat(1, facesList, idx);
    end

    S.Vertices = verticesList;
    S.Faces = facesList;
    S.FaceColor = 'none';
    S.EdgeColor = edgeColor;
    S.FaceAlpha = 0;
    S.LineWidth = lineWidth;
    patch(S, 'Parent', axesHandle);
end
