function generateNominalSpatialSamplingGrid(obj, samplingScheme)

    centerDegs = obj.theRGCMosaic.eccentricityDegs;
    sizeDegs = obj.theRGCMosaic.sizeDegs;
    gridHalfSamplesNum = 1;

    assert(numel(sizeDegs) == 2, 'MosaicPoolingOptimizer.generateNominalSpatialSamplingGrid: sizeDegs must be a 2-element vector');
    assert(ismember(samplingScheme, {'rectangular', 'hexagonal'}), 'MosaicPoolingOptimizer.generateNominalSpatialSamplingGrid: samplingScheme must be either ''hexagonal'' or ''rectangular''.');

    rgcRFpositionsDegs = obj.theRGCMosaic.rgcRFpositionsDegs;

    xCoords = squeeze(rgcRFpositionsDegs(:,1));
    yCoords = squeeze(rgcRFpositionsDegs(:,2));
    minX = min(xCoords(:)); maxX = max(xCoords(:));
    minY = min(yCoords(:)); maxY = max(yCoords(:));
    midX = (minX + maxX)/2;
    midY = (minY + maxY)/2;


    switch (samplingScheme)
        case 'rectangular'
            for iDim = 1:numel(sizeDegs)
               
                if (numel(gridHalfSamplesNum)==2)
                    nSamples = gridHalfSamplesNum(iDim);
                else
                    nSamples = gridHalfSamplesNum(1);
                end
        
                if (centerDegs(iDim) == 0)
                    x = logspace(log10(0.2/nSamples), log10(sizeDegs(iDim)*0.5), nSamples+1);
                    x(1) = 0;
                    coords{iDim} = unique([-fliplr(x) x]);
                else
                    x = linspace(0, sizeDegs(iDim)*0.5, nSamples);
                    coords{iDim} = unique([-fliplr(x) 0 x]);
                end
            end
            [X,Y] = meshgrid(coords{1}, coords{2});
            X = X(:);
            Y = Y(:);

        case 'hexagonal'
            if (gridHalfSamplesNum == 0)
                X = 0;
                Y = 0;
            else
                if (centerDegs(1) == 0) && (centerDegs(2) == 0)
                    if (gridHalfSamplesNum == 1)
                        r = logspace(log10(0.1), log10(max(sizeDegs)*0.3), gridHalfSamplesNum+1);
                    else
                        r = logspace(log10(0.1), log10(max(sizeDegs)*0.4), gridHalfSamplesNum+1);
                    end
                    r(1) = 0;
                    angles = -30+(0:60:360);
                    xScale = 0.75;
                    yScale = 1.0;
                else
                    xScale = 1;
                    yScale = 1;
                    r = linspace(0, max(sizeDegs)*0.4, gridHalfSamplesNum);
                    angles = 0:60:360; 
                end

                X = []; Y = [];
                delta = 0;
                for iRadius = numel(r):-1:1
                    if (delta == 0)
                        delta = 30;
                    else
                        delta = 0;
                    end
                    X(size(X,1)+(1:numel(angles)),1) = xScale * r(iRadius) * cosd(angles+delta);
                    Y(size(Y,1)+(1:numel(angles)),1) = yScale * r(iRadius) * sind(angles+delta);
                end
    
                % Add the 4 corner points
                if (centerDegs(1) == 0) && (centerDegs(2) == 0)
                    % Mosaic centered at (0,0)
                    X(size(X,1)+1,1) = sizeDegs(1)/2;
                    Y(size(Y,1)+1,1) = sizeDegs(2)/2;

                    X(size(X,1)+1,1) = sizeDegs(1)/2;
                    Y(size(Y,1)+1,1) =-sizeDegs(2)/2;

                    X(size(X,1)+1,1) =-sizeDegs(1)/2;
                    Y(size(Y,1)+1,1) = sizeDegs(2)/2;

                    X(size(X,1)+1,1) =-sizeDegs(1)/2;
                    Y(size(Y,1)+1,1) =-sizeDegs(2)/2;

                    % And 2 points along the y=0 meridian
                    X(size(X,1)+1,1) = sizeDegs(1)/2;
                    Y(size(Y,1)+1,1) = 0;

                    X(size(X,1)+1,1) = -sizeDegs(1)/2;
                    Y(size(Y,1)+1,1) = 0;



                else
                    % Off-center mosaic
                    X(size(X,1)+1,1) = minX-midX;
                    Y(size(Y,1)+1,1) = minY-midY;
                    
                    X(size(X,1)+1,1) = minX-midX; 
                    Y(size(Y,1)+1,1) = maxY-midY;

                    X(size(X,1)+1,1) = maxX-midX;
                    Y(size(Y,1)+1,1) = minY-midY;

                    X(size(X,1)+1,1) = maxX-midX;
                    Y(size(Y,1)+1,1) = maxY-midY;

                    % And 3 points along the midY line
                    X(size(X,1)+1,1) = minX-midX;
                    Y(size(Y,1)+1,1) = midY-midY;

                    X(size(X,1)+1,1) = maxX-midX;
                    Y(size(Y,1)+1,1) = midY-midY;

                    X(size(X,1)+1,1) = midX-midX;
                    Y(size(Y,1)+1,1) = midY-midY;
                end
            end  
    end % Switch


    if (centerDegs(1) == 0) && (centerDegs(2) == 0)
        if (max(abs(X(:)))>0)
            X = X(:)/max(abs(X(:))) * sizeDegs(1)*0.5*0.95;
        end
        if (max(abs(Y(:)))>0)
            Y = Y(:)/max(abs(Y(:))) * sizeDegs(2)*0.5*0.95;
        end
    end

    if (centerDegs(1) == 0) && (centerDegs(2) == 0)
        gridCoords(:,1) = X;
        gridCoords(:,2) = Y;
    else
        gridCoords(:,1) = midX + X;
        gridCoords(:,2) = midY + Y;
    end

    obj.nominalSpatialSamplingGrid = unique(gridCoords, 'rows');
end