function deconvolutionModel = computeDeconvolutionModel(obj)
    
    subjectsToInclude = [1:10]; %[4 8 9];
    tabulatedEccentricities = -[0 0.25 0.5 1 1.5 2:25];
    
    % Use WatsonRGCModel to retrieve cone apertures along the nasal meridian
    % for the tabulated eccentricities
    w = WatsonRGCModel();
    coneRFSpacingsDegs  = w.coneRFSpacingAndDensityAlongMeridian(abs(tabulatedEccentricities), ...
            'nasal meridian','deg', 'deg^2', ...
            'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', false);
    % Cone aperture is 80% of the cone radius
    coneApertureRadii = 0.8 * coneRFSpacingsDegs/2;
    coneApertureRadiusAtZeroDegs = coneApertureRadii(1);
    
    for eccIndex = 1:numel(tabulatedEccentricities)
        % Eccentricity
        eccDegs = tabulatedEccentricities(eccIndex);
        
        % Load deconvolution file for this eccentricity
        dataFileName = fullfile(obj.psfDeconvolutionDir,...
            sprintf('cellEcc_%2.1f_cellRefractionError_%2.2fD_VisualGain.mat', eccDegs, 0));
        load(dataFileName, 'retinalPoolingRadii', 'visualRadius', 'visualGain', 'subjectIDs', 'quadrants');
        retinalPoolingRadiiOriginal = retinalPoolingRadii;
        
        % Compute the minimum possible RF center radius as the mean (over all subjects)
        % convolution of the cone aperture at 0 deg with the PSF at 0 deg 
        if (eccDegs == 0)
            % Find the 2 closest retinal pooling radius
            [dd,idx] = sort(abs( coneApertureRadiusAtZeroDegs-retinalPoolingRadii));
            % Compute the minimal visual radius as the weighed average of the
            % visual radii corresponding to the 2 closest retinal pooling radii
            d1 = dd(2)/(dd(1)+dd(2));
            d2 = dd(1)/(dd(1)+dd(2));
            visualRadii = d1 * squeeze(visualRadius(idx(1),:,:)) + d2 * squeeze(visualRadius(idx(2),:,:));
            minVisualRadiusDegs = median(median(visualRadii));
        end
        
        % Get visual radius/gain data for the quadrant of interest
        quandrantsToInclude = 'both';   % choose from {'horizontal', 'vertical', 'both'}
        visualRadius = quadrantData(visualRadius, quandrantsToInclude, quadrants);
        visualGain = quadrantData(visualGain, quandrantsToInclude, quadrants);
        
        % We will fit the relationship between visual and retinal radii
        % using a saturating function, but only for retinal radii >= cone aperture radii
        % So select these data here:
        idx = find(retinalPoolingRadii >= coneApertureRadii(eccIndex));
        retinalPoolingRadii = retinalPoolingRadii(1,idx);
        visualRadius = visualRadius(:,idx);
        visualGain = visualGain(:,idx);
        
        % Compute the median visual radius over subjects/ecc quadrants
        visualRadius = visualRadius(subjectsToInclude,:);
        medianVisualRadius = (median(visualRadius,1, 'omitnan'))';
        
        % Compute the median gain over all subjects/ecc quadrants
        visualGain = visualGain(subjectsToInclude,:);
        medianVisualGain = (median(visualGain,1, 'omitnan'))';
       
        % Fit the relationship: retinalRadius(visualRadius) = Model(visualRadius, params)
        [modelFunctionRadius, fittedParamsRadius(eccIndex,:)] = fitRadiusData(medianVisualRadius, retinalPoolingRadii);
        % Fit the relationship: visualGain(retinalRadius) = Model(retinalRadius, params)
        [modelFunctionGain, fittedParamsGain(eccIndex,:)] = fitGainData(medianVisualGain, retinalPoolingRadii);
    end % eccIndex
    
    deconvolutionModel = struct(...
        'retinalPoolingRadii', retinalPoolingRadiiOriginal, ...
        'tabulatedEccentricities', tabulatedEccentricities, ...
        'fittedParamsRadius', fittedParamsRadius, ...
        'fittedParamsGain', fittedParamsGain, ...
        'modelFunctionRadius', modelFunctionRadius, ...
        'modelFunctionGain', modelFunctionGain, ...
        'minVisualRadiusDegs', minVisualRadiusDegs ...
        );
end

function [modelFunction, fittedParams] = fitRadiusData(visualRadius, retinalRadius)
    % For better fitting extend the visual radius data to 1.5
    visual  = [visualRadius; [0.7 0.8 1 1.5]'];
    retinal = [retinalRadius [0.7 0.8 1 1.5] ]';
        
    modelFunction = @(p,x)(p(1) - (p(2)-p(1))*exp(-p(3)*x));
    initialParams = [5 10 0.2];
    [fittedParams, fittedParamsSE] = nonLinearFitData(visual, retinal, modelFunction, initialParams);
end

function [modelFunction, fittedParams] = fitGainData(visualGain, retinalRadius)
    % For better fitting extend the visual radius data to 1.5
    visual  = [visualGain; [1 1 1]'];
    retinal = [retinalRadius [0.8 1 2]]';
        
    modelFunction = @(p,x)((p(1)*x.^p(3))./(x.^p(3)+p(2)));
    initialParams = [0.9 0.002 2];
    [fittedParams, fittedParamsSE] = nonLinearFitData(retinal, visual, modelFunction, initialParams);
end

function [fittedParams, fittedParamsSE] = nonLinearFitData(x,y, modelFunction, initialParams)
    opts.RobustWgtFun = 'talwar';
    [fittedParams,~,~,varCovarianceMatrix,~] = nlinfit(x,y,modelFunction,initialParams,opts);
    % standard error of the mean
    fittedParamsSE = sqrt(diag(varCovarianceMatrix));
    fittedParamsSE = fittedParamsSE';
end

function data = quadrantData(allQuadrantData, quandrantsToInclude, quadrants)
     switch (quandrantsToInclude)
         case  'horizontal'
             % Use horizontal eccs only
             data = squeeze(allQuadrantData(:, 1, :));
         case 'vertical'
             % Use vertical eccs only
             data =  [squeeze(allQuadrantData(:, 2, :)) squeeze(allQuadrantData(:, 3, :))];
         case 'both'
             data = [squeeze(allQuadrantData(:, 1, :))  squeeze(allQuadrantData(:, 2, :)) squeeze(allQuadrantData(:, 3, :))];
         otherwise
             error('quandrantsToInclude has an invalid value')
     end
     data = data';
 end