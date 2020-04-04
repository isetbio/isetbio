function visualizeRFs(connectivityMatrix, conePositionsMicrons, RGCRFPositionsMicrons, coneSpacingsMicrons, roi)

    % Generate RFs of RGCs based on cone psotions and connection matrix
    [RGCrfs, RGCRFsupportX, RGCRFsupportY] = generateRGCRFsFromConnectivityMatrix(...
        connectivityMatrix, conePositionsMicrons, coneSpacingsMicrons, roi);
  
    zLevels = 0.25 : 0.25 : 0.9;
    whichLevelsToContour = 1:numel(zLevels);
    
    hFig = figure(1); clf;
    theAxesGrid = plotlab.axesGrid(hFig, 'leftMargin', 0.03);
    theAxesGrid = theAxesGrid{1,1};
    set(theAxesGrid, 'XLim', roi.center(1)+1.1*roi.size(1)/2*[-1 1], 'YLim', roi.center(2)+1.1*roi.size(2)/2*[-1 1]);
    hold(theAxesGrid, 'on');
    
    rgcsNum = size(RGCrfs,1);
    for RGCindex = 1:rgcsNum
        theRF = squeeze(RGCrfs(RGCindex,:,:));
        if (all(isnan(theRF(:))))
            fprintf('Zero RF. skipping\n');
            continue;
        end
        contourf(RGCRFsupportX, RGCRFsupportY,theRF, zLevels);
        %C = contourc(RGCRFsupportX, RGCRFsupportY,theRF, zLevels);
        %renderContourPlot(theAxesGrid, C, zLevels, whichLevelsToContour );
    end
    colormap(brewermap(512, 'greys'))
    set(theAxesGrid, 'CLim', [0 1]);
end

function renderContourPlot(theAxes, C, zLevels, whichLevelsToContour )
    k = 1;
    while k < size(C,2)
        level = C(1,k);
        points = C(2,k);
        for kLevel = 1:numel(zLevels)
            if (level == zLevels(kLevel))
                xRGCEnsembleOutline.level{kLevel} = C(1,k+(1:points));
                yRGCEnsembleOutline.level{kLevel} = C(2,k+(1:points));
            end
        end
        k = k+points+1;
    end

    maxLevel = max(whichLevelsToContour);
    minLevel = min(whichLevelsToContour);
    
    for iLevel = 1:numel(whichLevelsToContour)
        theLevel = whichLevelsToContour(iLevel);
        
        if (theLevel == minLevel)
            edgeAlpha = 1;
        else
            edgeAlpha = 0;
        end
        
        f = 1:numel(xRGCEnsembleOutline.level{theLevel});
        v = [xRGCEnsembleOutline.level{theLevel}(:) yRGCEnsembleOutline.level{theLevel}(:)];
        patch(theAxes, 'Faces', f, 'Vertices', v, 'FaceColor', (1-0.8*theLevel/maxLevel)*[0.8 0.8 0.8], ...
            'FaceAlpha', 1, 'EdgeColor', [0 0 0.6], 'EdgeAlpha',edgeAlpha, 'LineWidth', 1.0);

    end
end


function [RGCrfs, xAxis, yAxis] = generateRGCRFsFromConnectivityMatrix(connectivityMatrix, conePositions, coneSpacings,  roi)
    % Sampling for contours
    deltaX = 0.2;
    xAxis = (roi.center(1)-roi.size(1)/2): deltaX: (roi.center(1)+roi.size(1)/2);
    yAxis = (roi.center(2)-roi.size(2)/2): deltaX: (roi.center(2)+roi.size(2)/2);
    [X,Y] = meshgrid(xAxis,yAxis);
    
    conesNum = size(connectivityMatrix,1);
    rgcsNum = size(connectivityMatrix,2);
    
    coneProfiles = zeros(conesNum, size(X,1), size(X,2));
    
    for coneIndex = 1:conesNum 
        cP = squeeze(conePositions(coneIndex,:));
        coneSigma = coneSpacings(coneIndex)/6;
        
        coneProfile = exp(-0.5*((X-cP(1))/coneSigma).^2) .* exp(-0.5*((Y-cP(2))/coneSigma).^2);
        coneProfile = coneProfile.^0.2;
        coneProfiles(coneIndex,:,:) = coneProfile/max(coneProfile(:));
    end
    
    RGCrfs = zeros(rgcsNum, size(X,1), size(X,2));
   
    for RGCindex = 1:rgcsNum
        RGCrf = squeeze(RGCrfs(RGCindex,:,:));
        for coneIndex = 1:conesNum
            if (connectivityMatrix(coneIndex, RGCindex)>0)
                RGCrf = RGCrf + squeeze(coneProfiles(coneIndex,:,:)) * connectivityMatrix(coneIndex, RGCindex);
            end
        end
        RGCrfs(RGCindex,:,:) = RGCrf / max(RGCrf(:));
    end
end




