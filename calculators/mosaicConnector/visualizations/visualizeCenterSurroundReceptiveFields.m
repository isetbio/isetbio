function visualizeCenterSurroundReceptiveFields(figNo, synthesizedRFParams, dataSet, ...
    conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
    RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
    midgetRGCconnectionMatrix, ...
    plotlabOBJ, pdfFileName, exportsDir)

    % Get data
    switch (dataSet)
        case 'visual'
            rfParams = synthesizedRFParams.visual;
        case 'retinal'
            rfParams = synthesizedRFParams.retinal;
        otherwise
            error('Unknown dataSet: ''%s''.', dataSet)
    end
    
    % Ecc and indices of RGCs in the analyzed patch
    eccDegs = synthesizedRFParams.eccDegs;
    rgcIndices = synthesizedRFParams.rgcIndices;
    
    % Analyzed patch
    patchEccDegs = synthesizedRFParams.patchEccDegs;
    patchSizeDegs = synthesizedRFParams.patchSizeDegs;

    hFig = figure();
    
    xx = cosd(0:5:360);
    yy = sind(0:5:360);
    for iRGC = 1:numel(rgcIndices)

        rgcIndex = synthesizedRFParams.rgcIndices(iRGC);
        rgcPositionMicrons = synthesizedRFParams.rgcCenterPositionMicrons(iRGC,:);
        rgcCenterRadiusDegs = rfParams.centerRadiiDegs(iRGC);
        rgcSurroundRadiusDegs = rfParams.surroundRadiiDegs(iRGC);
        rgcCenterRadiusMicrons = convertSizeDegsToSizeMicrons(rgcCenterRadiusDegs, eccDegs(iRGC));
        rgcSurroundRadiusMicrons = convertSizeDegsToSizeMicrons(rgcSurroundRadiusDegs, eccDegs(iRGC));
       
        % [HERE]
        % Determine cone weights to the surround region
        
        % Center connectivity vector
        connectivityVector = full(squeeze(midgetRGCconnectionMatrix(:, rgcIndex)));
        connectedConeIDs = find(connectivityVector>0);
        
        % Plot cones connected to the RF center
        for iCone = 1:numel(connectedConeIDs)
            coneIndex = connectedConeIDs(iCone);
            switch (coneTypes(coneIndex))
                case 2
                    coneColor = [1 0 0];
                case 3
                    coneColor = [0 1 0];
                case 4
                    coneColor = [0 0 1];
            end
            cPosMicrons = conePositionsMicrons(coneIndex,:);
            cRadiusMicrons = 0.8*0.5*coneSpacingsMicrons(coneIndex);
            cProfile(:,1) = cPosMicrons(1) + cRadiusMicrons*xx;
            cProfile(:,2) = cPosMicrons(2) + cRadiusMicrons*yy;
            patch('faces', 1:size(cProfile,1), 'Vertices', cProfile, ...
                'FaceColor', 0.5*(coneColor+[0.5 0.5 0.5]), 'FaceAlpha', 1, 'EdgeColor', coneColor, 'EdgeAlpha', 0.8);
        
            hold on;
        end
        
        % Plot the outline of the surround (at 1/e)
        rgcSurroundProfile(:,1) = rgcPositionMicrons(1) + rgcSurroundRadiusMicrons * xx;
        rgcSurroundProfile(:,2) = rgcPositionMicrons(2) + rgcSurroundRadiusMicrons * yy;
        patch('faces', 1:size(rgcSurroundProfile,1), 'Vertices', rgcSurroundProfile, ...
            'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeColor', [0 0 0], 'EdgeAlpha', 0.8);
        
        % Plot the outline of the center (at 1/e)
        rgcCenterProfile(:,1) = rgcPositionMicrons(1) + rgcCenterRadiusMicrons * xx;
        rgcCenterProfile(:,2) = rgcPositionMicrons(2) + rgcCenterRadiusMicrons * yy;
        plot(rgcCenterProfile(:,1), rgcCenterProfile(:,2), 'k-');
        
        
        
        axis 'equal'
        drawnow;
        pause
    end % iRGC
    
    plotlabOBJ.exportFig(hFig, 'png', sprintf('%s__%s',pdfFileName, dataSet),exportsDir);
end

function sizeMicrons = convertSizeDegsToSizeMicrons(sizeDegs, eccDegs)
    sizeMicrons = ...
        1000*(WatsonRGCModel.rhoDegsToMMs(eccDegs + sizeDegs/2) - WatsonRGCModel.rhoDegsToMMs(eccDegs - sizeDegs/2));
end
