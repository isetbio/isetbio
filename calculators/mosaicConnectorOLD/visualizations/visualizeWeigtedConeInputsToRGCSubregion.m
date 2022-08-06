function visualizeWeigtedConeInputsToRGCSubregion(theAxes, coneTypes, ...
    conePositionsMicrons, coneSpacingsMicrons, connectionWeights,  maxWeightVisualized, ...
    rgcPositionMicrons,  spatialSupportRangeMicrons, polarity)

    xx = cosd(0:5:360);
    yy = sind(0:5:360);
    
    realConnectionWeights = connectionWeights;
    
    % For visualization we scale these down
    connectionWeights = connectionWeights / maxWeightVisualized;
    connectionWeights = min([ones(numel(connectionWeights),1) connectionWeights],[],2);
    
    backgroundColor = [1 1 1];
    cla(theAxes);
    
    x = -spatialSupportRangeMicrons:0.5:spatialSupportRangeMicrons;
    coneSigma = 0.8*0.5*mean(coneSpacingsMicrons)/3;
    coneKernel = (exp(-0.5*(x/coneSigma).^2)).^0.5;
    
    y = x * 0;
    for iCone = 1:numel(connectionWeights)
        switch (coneTypes(iCone))
            case 2
               coneColor = [1 0 0]*connectionWeights(iCone) + backgroundColor * (1-connectionWeights(iCone));
            case 3
               coneColor = [0 0.7 0]*connectionWeights(iCone) + backgroundColor * (1-connectionWeights(iCone));
            case 4
               coneColor = [0 0 1]*connectionWeights(iCone) + backgroundColor * (1-connectionWeights(iCone));
        end
        coneColor = min([ones(1,3); coneColor]);
        cPosMicrons = conePositionsMicrons(iCone,:);
        cRadiusMicrons = 0.8*0.5*coneSpacingsMicrons(iCone);
        cProfile(:,1) = cPosMicrons(1) + cRadiusMicrons*xx;
        cProfile(:,2) = cPosMicrons(2) + cRadiusMicrons*yy;
        xCoord = cPosMicrons(1)-rgcPositionMicrons(1)+spatialSupportRangeMicrons;
        if (xCoord>0)&&(xCoord<=2*spatialSupportRangeMicrons)
            y(round(2*xCoord)+1) = y(round(2*xCoord)+1) + realConnectionWeights(iCone);
        end
        patch(theAxes, 'faces', 1:size(cProfile,1), 'Vertices', cProfile, ...
                'FaceColor', coneColor, 'FaceAlpha', 1, 'EdgeColor', [0.4 0.4 0.4], 'EdgeAlpha', 0.4);
        hold(theAxes, 'on');
    end % iCone
    
    visualizationGain = 0.05;
    y = visualizationGain*conv(y, coneKernel, 'same');
    baseline = -2*spatialSupportRangeMicrons/3+rgcPositionMicrons(2);
    y = baseline + y*polarity;
    area(theAxes, rgcPositionMicrons(1)+x, y, baseline, 'FaceColor',[ 0.5 0.5 0.5] + polarity*[0.3 0.3 0.3], 'EdgeAlpha', 0.2);
    
    % Plot the outline of the subregion
%     plotSubregionOutline = true;
%     if (plotSubregionOutline)
%         rgcSubregionOutline(:,1) = rgcPositionMicrons(1) + rgcSubregionRadiusMicrons * xx;
%         rgcSubregionOutline(:,2) = rgcPositionMicrons(2) + rgcSubregionRadiusMicrons * yy;
% 
%         patch(theAxes,'faces', 1:size(rgcSubregionOutline,1), 'Vertices', rgcSubregionOutline, ...
%                 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.0, 'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.3);
%     end
    grid(theAxes, 'off');
    set(theAxes, 'XLim', rgcPositionMicrons(1) + spatialSupportRangeMicrons*[-1 1], ...
                 'YLim', rgcPositionMicrons(2) + spatialSupportRangeMicrons*[-1 1], ...
                 'XTick', round(rgcPositionMicrons(1)), 'YTick', round(rgcPositionMicrons(2)), ...
                 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
    
    axis(theAxes,'square');     
end