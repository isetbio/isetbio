function gridCoords = eccentricitySamplingGridCoords(eccentricityDegs, sizeDegs, gridHalfSamplesNum, samplingScheme, visualizeGridCoords)

    assert(numel(sizeDegs) == 2, 'MidgetRGCMosaic:gridNormalizedCoords: sizeDegs must be a 2-element vector');
    assert(ismember(samplingScheme, {'rectangular', 'hexagonal'}), 'MidgetRGCMosaic:gridNormalizedCoords: samplingScheme must be either ''hexagonal'' or ''rectangular''.');


    switch (samplingScheme)
        case 'rectangular'
            for iDim = 1:numel(sizeDegs)
               
                if (numel(gridHalfSamplesNum)==2)
                    nSamples = gridHalfSamplesNum(iDim);
                else
                    nSamples = gridHalfSamplesNum(1);
                end
        
                if (eccentricityDegs(iDim) == 0)
                    x = logspace(log10(0.2/nSamples), log10(sizeDegs(iDim)*0.5), nSamples+1);
                    x(1) = 0;
                    coords{iDim} = unique([-fliplr(x) x]);
                else
                    x = linspace(0, sizeDegs(iDim)*0.5, nSamples);
                    coords{iDim} = unique([-fliplr(x) 0 x]);
                end
            end
            [X,Y] = meshgrid(coords{1}, coords{2});
    
        case 'hexagonal'
            if (gridHalfSamplesNum == 0)
                X = 0;
                Y = 0;
            else
                r = linspace(0, max(sizeDegs)*0.5, gridHalfSamplesNum+1);
                angles = 0:60:360;
                X = [];
                Y = [];
        
                delta = 0;
                for iRadius = numel(r):-1:1
                    if (delta == 0)
                        delta = 30;
                    else
                        delta = 0;
                    end
                    X(size(X,1)+(1:numel(angles)),1) = r(iRadius) * cosd(angles+delta);
                    Y(size(Y,1)+(1:numel(angles)),1) = r(iRadius) * sind(angles+delta);
                end
    
            end
    end


    if (max(abs(X(:)))>0)
        X = X(:)/max(abs(X(:))) * sizeDegs(1)*0.5;
    end
    if (max(abs(Y(:)))>0)
        Y = Y(:)/max(abs(Y(:))) * sizeDegs(2)*0.5;
    end

    gridCoords(:,1) = eccentricityDegs(1) + X;
    gridCoords(:,2) = eccentricityDegs(2) + Y;
    gridCoords = unique(gridCoords, 'rows');

    if (visualizeGridCoords)
        hFig = figure(1); clf;
        set(hFig, 'Color', [1 1 1]);
        
        xx = eccentricityDegs(1) + sizeDegs(1)/2*[-1 -1 1  1 -1];
        yy = eccentricityDegs(2) + sizeDegs(2)/2*[-1  1 1 -1 -1];
        plot(xx,yy,'k-', 'LineWidth', 1);
        hold on;
        plot(gridCoords(:,1), gridCoords(:,2), 'r+', 'MarkerSize', 12, 'LineWidth', 3.0);
        set(gca, 'FontSize', 16);
        axis 'equal'
        box off;
        grid on
        set(gca, 'XLim', eccentricityDegs(1) + sizeDegs(1)*0.5*[-1 1] + [-0.1 0.1], 'YLim', eccentricityDegs(2) + sizeDegs(2)*0.5*[-1 1] + [-0.1 0.1]);
    end
end
