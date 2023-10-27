function contourData = subregionEllipseFromPooledCones(...
    conePos, coneRc, poolingWeights, ...
    xSupport, ySupport, spatialSupportSamples, centerSubregionContourSamples)

    % Compute spatial support
    xSep = max(coneRc)*2.5*sqrt(numel(poolingWeights));

    if (isempty(poolingWeights))
        fprintf(2, 'poolingWeights is []. Returning an empty contourData struct\n');
        contourData = [];
        return;
    end

    if (isempty(xSupport))
        xx = conePos(:,1);
        xSupport = linspace(min(xx)-xSep,max(xx)+xSep,spatialSupportSamples);
    end

    if (isempty(ySupport))
        yy = conePos(:,2);
        ySupport = linspace(min(yy)-xSep,max(yy)+xSep,spatialSupportSamples);
    end

    [X,Y] = meshgrid(xSupport, ySupport);
    RF = zeros(size(X));

    if (numel(poolingWeights)>1)
        xx = conePos(:,1);
        yy = conePos(:,2);
        dx = max(xx)-min(xx);
        dy = max(yy)-min(yy);
        if (dx > dy)
            [~,idx] = sort(xx);
        else
            [~,idx] = sort(yy);
        end
        
        xo = mean(xx);
        yo = mean(yy);

        poolingWeights = poolingWeights(idx);
        coneRc = coneRc(idx);
        xx = xx(idx);
        yy = yy(idx);
        interpolationPointsNum = 4;
        xxxInterp = zeros(1, interpolationPointsNum * numel(xx));
        yyyInterp = zeros(1, interpolationPointsNum * numel(xx));
        coneRcInterp = zeros(1, interpolationPointsNum * numel(xx));
        poolingWeightsInterp = zeros(1, interpolationPointsNum * numel(xx));

        for ii = 1:numel(xx)
            for iii = 1:interpolationPointsNum
                f = (iii-1)/interpolationPointsNum;
                xxxInterp((ii-1)*interpolationPointsNum+iii) = xo*(1-f) + xx(ii)*f;
                yyyInterp((ii-1)*interpolationPointsNum+iii) = yo*(1-f) + yy(ii)*f;
                coneRcInterp((ii-1)*interpolationPointsNum+iii) = coneRc(ii);
                poolingWeightsInterp((ii-1)*interpolationPointsNum+iii) = poolingWeights(ii);
            end
        end


        for iConeInterp = 1:numel(coneRcInterp)
            % Characteristic radius of the input RF
            rC = coneRcInterp(iConeInterp);
            % Compute aperture2D x weight
            XX = X-xxxInterp(iConeInterp);
            YY = Y-yyyInterp(iConeInterp);
            theAperture2D = poolingWeightsInterp(iConeInterp) * exp(-(XX/rC).^2) .* exp(-(YY/rC).^2);
            % Accumulate 2D apertures
            RF = RF + theAperture2D;
            RF(RF>1) = 1;
        end

        % for iCone = 1:numel(poolingWeights)
        %     % Characteristic radius of the input RF
        %     rC = coneRc(iCone);
        %     % Compute aperture2D x weight
        %     XX = X-conePos(iCone,1);
        %     YY = Y-conePos(iCone,2);
        %     theAperture2D = poolingWeights(iCone) * exp(-(XX/rC).^2) .* exp(-(YY/rC).^2);
        %     % Accumulate 2D apertures
        %     RF = RF + theAperture2D;
        % end
    else
        % Characteristic radius of the input RF
        iCone = 1;
        rC = coneRc(iCone);
        % Compute aperture2D x weight
        XX = X-conePos(iCone,1);
        YY = Y-conePos(iCone,2);
        theAperture2D = poolingWeights(iCone) * exp(-(XX/rC).^2) .* exp(-(YY/rC).^2);
        RF = theAperture2D;
        RF(RF>1) = 1;
    end


    zLevel = 0.02;
    contourData = mRGCMosaic.ellipseContourFromSubregionRFmap(xSupport, ySupport, RF, zLevel, centerSubregionContourSamples);


    % % Binarize
    % RF = RF / max(RF(:));
    % RF(RF<0.1) = 0.0;
    % RF(RF>0) = 1.0;
    % BW = imbinarize(RF);
    % 
    % 
    % % Extract the maximum area
    % BW = imclearborder(BW);
    % BW = bwareafilt(BW,1);
    % 
    % % Calculate centroid, orientation and major/minor axis length of the ellipse
    % s = regionprops(BW,{'Centroid','Orientation','MajorAxisLength','MinorAxisLength'});
    % if (isempty(s))
    %     figure()
    %     subplot(1,2,1);
    %     imagesc(RF)
    %     axis 'image'
    %     subplot(1,2,2)
    %     imagesc(BW);
    %     axis 'image'
    % 
    %     contourData = [];
    %     return;
    % end
    % 
    % % Calculate the ellipse line
    % theta = linspace(0, 2*pi, centerSubregionContourSamples);
    % col = (s.MajorAxisLength/2)*cos(theta);
    % row = (s.MinorAxisLength/2)*sin(theta);
    % M = makehgtform('translate',[s.Centroid, 0],'zrotate',deg2rad(-1*s.Orientation));
    % D = M*[col;row;zeros(1,numel(row));ones(1,numel(row))];
    % 
    % x = D(1,:);
    % y = D(2,:);
    % x = (x-1)/(numel(xSupport)-1) * (xSupport(end)-xSupport(1)) + xSupport(1); 
    % y = (y-1)/(numel(ySupport)-1) * (ySupport(end)-ySupport(1)) + ySupport(1); 
    % 
    % v = [x(:) y(:)];
    % f = 1:numel(x);
    % s = struct('faces', f, 'vertices', v);
    % contourData = s;
end