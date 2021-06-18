classdef constants
% Constants from Polans study

    properties (Constant)
        % Label for left eye
        leftEye = 'left eye';
        
        % Label for right eye
        rightEye = 'right eye';
        
        % Label for nasal meridian
        nasalMeridian    = 'nasal meridian';
        
        % Label for temporal meridian
        temporalMeridian = 'temporal meridian';
        
        % Label for inferior meridian
        inferiorMeridian = 'inferior meridian';
        
        % Label for superior meridian
        superiorMeridian = 'superior meridian';
        
        % Pupil diameter for all measurements
        measurementPupilDiamMM = 4.0;
        
        % All data are for the right eye
        measurementEye = 'right';
        
        % Eccentricities at which measurements were made       
        measurementHorizontalEccentricities = 40:-1:-40;
        measurementVerticalEccentricities = 25:-5:-25;
    end
   
    methods (Static)
        function retinalQuadrantName = retinalQuadrantFromEcc(ecc, whichEye)
            assert(numel(ecc) == 2, 'Eccentricity must be a 2-element vector (x,y)');
            assert(ismember(whichEye, {PolansOptics.constants.leftEye, PolansOptics.constants.leftEye}));
            switch (whichEye)
                case 'right'
                    if (ecc(1) < 0)
                        horizontalMeridianName = PolansOptics.constants.nasalMeridian;
                    elseif (ecc(1) == 0)
                        horizontalMeridianName = 'foveal pit';
                    else
                        horizontalMeridianName = PolansOptics.constants.temporalMeridian;
                    end
                case 'left'
                    if (ecc(1) < 0)
                        horizontalMeridianName = PolansOptics.constants.temporalMeridian;
                    elseif (ecc(1) == 0)
                        horizontalMeridianName = 'foveal pit';
                    else
                        horizontalMeridianName = PolansOptics.constants.nasalMeridian;
                    end
            end
            if (ecc(2) < 0)
                verticalMeridianName = PolansOptics.constants.inferiorMeridian;
            elseif (ecc(1) == 0)
                verticalMeridianName = 'foveal pit';
            else
                verticalMeridianName = PolansOptics.constants.superiorMeridian;
            end
            retinalQuadrantName{1} = horizontalMeridianName;
            retinalQuadrantName{2} = verticalMeridianName;
        end
        
        function [zMap, zCoeffIndices] = ZernikeCoefficientsMap(subjectIndex)
            % Import raw data
            allData = rawDataReadData('zCoefsPolans2015', 'datatype', 'isetbiomatfileonpath');
            allData = allData.data;
    
            % Reported Z-coeffs are Z3-Z20
            zCoeffIndices = 3:size(allData,3);
            zCoeffs = squeeze(allData(subjectIndex , :, zCoeffIndices));
            
            % Assemble xy map of Zernike coefficients
            zMap = zeros(...
                numel(PolansOptics.constants.measurementVerticalEccentricities), ...
                numel(PolansOptics.constants.measurementHorizontalEccentricities), ...
                numel(zCoeffIndices));
            
            xyIndex = 0;
            for vEccIndex = 1:size(zMap,1)
                for hEccIndex = 1:size(zMap,2)
                    xyIndex = xyIndex+1;
                    zMap(vEccIndex, hEccIndex,:) = zCoeffs(xyIndex,:);
                end
            end
        end
        
        function [zMapDeNoised, zCoeffIndices, outlierPointsMap] = deNoisedZernikeCoefficientsMap(subjectIndex, movingMeanWindowSize, sigmaMultiplier)
            % Obtain raw data set
            [zMap, zCoeffIndices] = PolansOptics.constants.ZernikeCoefficientsMap(subjectIndex);
            zMapDeNoised = zMap*0;
            outlierPointsMap = zMap*0;
            
            fprintf('DeNoised Zcoeffs: moving window size = %d, sigmaMultipiler = %2.2f\n', movingMeanWindowSize, sigmaMultiplier);
            for eccYindex = 1:numel(PolansOptics.constants.measurementVerticalEccentricities)
            for eccXindex = 1:numel(PolansOptics.constants.measurementHorizontalEccentricities)
                zCoeffs = squeeze(zMap(eccYindex,eccXindex,:));
                if (all(zCoeffs == 0))
                    zMap(eccYindex,eccXindex,:) = nan;
                end
            end
            end
        
            % Fit polynomials to obtain smooth Zernike coefficients for all eccentricities
            for zCoeffIndex = 1:size(zMap,3)
                for eccYindex = 1:numel(PolansOptics.constants.measurementVerticalEccentricities)

                    % Raw z-coeffs
                    zCoeffs = squeeze(zMap(eccYindex,:, zCoeffIndex));
                    
                    % Find the trend of z-coeffs using a moving average window
                    trend = movmean(zCoeffs,movingMeanWindowSize, 'omitnan');
                    
                    % Remove the trend to find outliers
                    detrendedZCoeffs = zCoeffs - trend;
                    
                    % Indices of outliers (far away from the trend
                    % (sigmaMultiplier * std) or having nan
                    idx = (abs(detrendedZCoeffs) > sigmaMultiplier * std(detrendedZCoeffs, 'omitnan')) | isnan(zCoeffs);
                    
                    % Keep track of outlier points (useful for plotting them)
                    outlierPointsMap(eccYindex, :, zCoeffIndex) = idx;

                    % Indices to use for polynomial fit (all except the outliers)
                    indicesForPolynomialFit = find(idx == 0);

                    % Source data (excluding noisy data points)
                    rawDataX = PolansOptics.constants.measurementHorizontalEccentricities(indicesForPolynomialFit);
                    rawDataY = zMap(eccYindex,indicesForPolynomialFit, zCoeffIndex);
                    
                    % Find bet fit polynomial up
                    maxOrder = 4;
                    if (zCoeffIndex > 3)
                        maxOrder = 6;
                    end
                    examinedOrders = 2:maxOrder;
                    residuals = inf(1,maxOrder);
                    p = cell(1,  maxOrder);

                    for k = 1:numel(examinedOrders)
                        theOrder = examinedOrders(k);
                        p{k} = polyfit(rawDataX, rawDataY,theOrder); 
                        defocusCoeffsMicronsInterpolated = polyval(p{k} ,PolansOptics.constants.measurementHorizontalEccentricities);
                        differences = zMap(eccYindex,indicesForPolynomialFit, zCoeffIndex) - defocusCoeffsMicronsInterpolated(indicesForPolynomialFit);
                        residuals(k) = sum(differences.^2);
                    end
                    
                    % Find the minimum error and the corresponding polynomial order
                    [~,theBestK] = min(residuals);
                    
                    % Polynomial fit for all eccentricities
                    intepolatedZoeffs = polyval(p{theBestK}, PolansOptics.constants.measurementHorizontalEccentricities);
                    
                    
                    figure(100); clf;
                    plot(PolansOptics.constants.measurementHorizontalEccentricities, zCoeffs, 'ks', 'MarkerSize', 12);
                    hold on;
                    plot(PolansOptics.constants.measurementHorizontalEccentricities, intepolatedZoeffs, 'r-', 'LineWidth', 1.5);
                    drawnow;
                    title(sprintf('Z%d', zCoeffIndex));
                    pause(1.0)
                    
                    doSecondPass = false;
                    if (doSecondPass) 
                        % Second pass
                        residuals = zCoeffs - intepolatedZoeffs;
                        indicesForPolynomialFit = find((abs(residuals) < sigmaMultiplier * std(residuals, 'omitnan')) &  (~isnan(zCoeffs)));


                        % Source data (excluding noisy data points)
                        rawDataX = PolansOptics.constants.measurementHorizontalEccentricities(indicesForPolynomialFit);
                        rawDataY = zMap(eccYindex,indicesForPolynomialFit, zCoeffIndex);

                        % Find best fit polynomial up to order maxOrder
                        examinedOrders = 2:maxOrder;
                        residuals = inf(1,maxOrder);
                        p = cell(1,  maxOrder);

                        for k = 1:numel(examinedOrders)
                            theOrder = examinedOrders(k);
                            p{k} = polyfit(rawDataX, rawDataY,theOrder); 
                            defocusCoeffsMicronsInterpolated = polyval(p{k} ,PolansOptics.constants.measurementHorizontalEccentricities);
                            differences = zMap(eccYindex,indicesForPolynomialFit, zCoeffIndex) - defocusCoeffsMicronsInterpolated(indicesForPolynomialFit);
                            residuals(k) = sum(differences.^2);
                        end

                        % Find the minimum error, and the corresponding
                        % polynomial order
                        [~,theBestK]= min(residuals);

                        % Polynomial fit for all eccentricities
                        intepolatedZoeffs = polyval(p{theBestK}, PolansOptics.constants.measurementHorizontalEccentricities);
                        end

                        % Denoised coeffs
                        zMapDeNoised(eccYindex,:, zCoeffIndex) = intepolatedZoeffs;
                end % eccYindex
            end % zCoeffIndex
        end
        
        
        
    end % Static methods
    
end
