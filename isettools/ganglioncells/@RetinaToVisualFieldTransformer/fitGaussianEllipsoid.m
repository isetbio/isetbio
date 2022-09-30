% Method to fit a 2D Gaussian ellipsoid
function theFittedGaussian = fitGaussianEllipsoid(supportX, supportY, theRF, varargin)

    p = inputParser;
    p.addParameter('flatTopGaussian', false, @islogical);
    p.addParameter('forcedOrientationDegs', [], @(x)(isempty(x) || isscalar(x)));
    p.addParameter('forcedEllipseRcYRcXratio', [], @(x)(isempty(x) || isscalar(x)));
    p.addParameter('rangeForEllipseRcYRcXratio', [], @(x)(isempty(x) || (isnumeric(x)&&(numel(x)==2)) ));
    p.addParameter('forcedCentroidXYpos', [], @(x)(isempty(x) || (isnumeric(x)&&(numel(x)==2)) ));
    p.addParameter('globalSearch', false, @islogical);
    p.addParameter('multiStartsNum', [], @(x)(isempty(x) || isscalar(x)));
    p.parse(varargin{:});
    flatTopGaussian = p.Results.flatTopGaussian;
    forcedEllipseRcYRcXratio = p.Results.forcedEllipseRcYRcXratio;
    forcedOrientationDegs = p.Results.forcedOrientationDegs;
    rangeForEllipseRcYRcXratio = p.Results.rangeForEllipseRcYRcXratio;
    forcedCentroidXYpos = p.Results.forcedCentroidXYpos;
    globalSearch = p.Results.globalSearch;
    multiStartsNum = p.Results.multiStartsNum;
    
    [X,Y] = meshgrid(supportX, supportY);
    xydata(:,:,1) = X;
    xydata(:,:,2) = Y;

    maxRF = max(theRF(:));
    theRF = theRF / maxRF;
    [theCentroid, theAxesLengths, theRotationAngle] = RetinaToVisualFieldTransformer.estimateGeometry(supportX, supportY, theRF);

    if (~isempty(forcedOrientationDegs))
        theRotationAngle = forcedOrientationDegs;
    end

    % Form parameter vector: [gain, xo, RcX, yo, RcY/RcX, rotationAngleDegs]
    params.initialValues = [...
        max(theRF(:)), ...
        theCentroid(1), ...
        theAxesLengths(1)/3, ...
        theCentroid(2), ...
        theAxesLengths(2)/theAxesLengths(1), ...
        theRotationAngle];

    % Form lower and upper value vectors
    params.lowerBounds = [ ...
        0 ...
        min(supportX) ...
        0*(max(supportX)-min(supportX)) ...
        min(supportY) ...
        1e-2  ...
        theRotationAngle-90];

    params.upperBounds = [ ...
        1 ...
        max(supportX) ...
        max(supportX)-min(supportX) ...
        max(supportY) ...
        1e2 ...
        theRotationAngle+90];

    if (~isempty(forcedEllipseRcYRcXratio))
        params.initialValues(5) = forcedEllipseRcYRcXratio;
        params.lowerBounds(5) = forcedEllipseRcYRcXratio;
        params.upperBounds(5) = forcedEllipseRcYRcXratio;
    else
        if (~isempty(rangeForEllipseRcYRcXratio))
            params.lowerBounds(5) = rangeForEllipseRcYRcXratio(1);
            params.upperBounds(5) = rangeForEllipseRcYRcXratio(2);
        end

    end

    if (~isempty(forcedCentroidXYpos))
        params.initialValues(2) = forcedCentroidXYpos(1);
        params.lowerBounds(2) = params.initialValues(2);
        params.upperBounds(2) = params.lowerBounds(2);

        params.initialValues(4) = forcedCentroidXYpos(2);
        params.lowerBounds(4) = params.initialValues(4);
        params.upperBounds(4) = params.lowerBounds(4);
    end

    if (~isempty(forcedOrientationDegs))
        params.lowerBounds(end) = forcedOrientationDegs;
        params.upperBounds(end) = forcedOrientationDegs;
    end

    if (flatTopGaussian)
        params.initialValues(numel(params.initialValues)+1) = 1.5;
        params.lowerBounds(numel(params.lowerBounds)+1) = 1;
        params.upperBounds(numel(params.upperBounds)+1) = 2;

        params.initialValues(numel(params.initialValues)+1) = 1.5;
        params.lowerBounds(numel(params.lowerBounds)+1) = 1;
        params.upperBounds(numel(params.upperBounds)+1) = 2;

    else
        params.initialValues(numel(params.initialValues)+1) = 1.0;
        params.lowerBounds(numel(params.lowerBounds)+1) = 1.0;
        params.upperBounds(numel(params.upperBounds)+1) = 1.0;

        params.initialValues(numel(params.initialValues)+1) = 1.0;
        params.lowerBounds(numel(params.lowerBounds)+1) = 1.0;
        params.upperBounds(numel(params.upperBounds)+1) = 1.0;

    end

    % Do the fitting
    if (globalSearch)

        % Ready to fit
        options = optimset(...
            'Display', 'off', ...
            'Algorithm', 'interior-point',... % 'sqp', ... % 'interior-point',...
            'GradObj', 'off', ...
            'DerivativeCheck', 'off', ...
            'MaxFunEvals', 10^5, ...
            'MaxIter', 10^3);

        % Multi-start
        problem = createOptimProblem('fmincon',...
              'objective', @gaussian2DObjective, ...
              'x0', params.initialValues, ...
              'lb', params.lowerBounds, ...
              'ub', params.upperBounds, ...
              'options', options...
          );

         ms = MultiStart(...
              'Display', 'off', ...
              'StartPointsToRun','bounds-ineqs', ...  % run only initial points that are feasible with respect to bounds and inequality constraints.
              'UseParallel', true);
      
         % Run the multi-start
         if (isempty(multiStartsNum))
            multiStartsNum = 32;
         end

         fittedParams = run(ms, problem, multiStartsNum);
    else
        % Local search
        [fittedParams,resnorm,residual,exitflag] = lsqcurvefit(@gaussian2D,params.initialValues,xydata,theRF,params.lowerBounds,params.upperBounds);
    end




    xo = fittedParams(2);
    yo = fittedParams(4);
    RcX = fittedParams(3);
    RcY = RcX * fittedParams(5);
 

    % Compute the fitted 2D Gaussian
    theFittedGaussian.ellipsoidMap = gaussian2D(fittedParams,xydata);
    theFittedGaussian.characteristicRadii = [RcX RcY];
    theFittedGaussian.rotationDegs = fittedParams(6);
    theFittedGaussian.flatTopExponents = [fittedParams(7) fittedParams(8)];
    theFittedGaussian.xyCenter = [xo yo];
    

    % Nested function gaussian2DObjective
    function rmse = gaussian2DObjective(params)
        fittedRF = gaussian2D(params, xydata);
        fullRMSE = ((fittedRF(:) - theRF(:))).^2;
        rmse =  sqrt(mean(fullRMSE,1));
    end

end

function F = gaussian2D(params,xydata)
    % Retrieve spatial support
    X = squeeze(xydata(:,:,1));
    Y = squeeze(xydata(:,:,2));

    % Retrieve params
    gain = params(1);
    xo = params(2);
    yo = params(4);
    RcX = params(3);
    RcY = RcX*params(5);
    rotationAngle = params(6);
    flatTopGaussianExponentX = params(7);
    flatTopGaussianExponentY = params(8);

    % Apply axes rotation
    Xrot = X * cosd(rotationAngle) -  Y*sind(rotationAngle);
    Yrot = X * sind(rotationAngle) +  Y*cosd(rotationAngle);
    xorot = xo * cosd(rotationAngle) -  yo*sind(rotationAngle);
    yorot = xo * sind(rotationAngle) +  yo*cosd(rotationAngle);

    % Compute 2D Gaussian
    exponentX = 2.0 * flatTopGaussianExponentX;
    exponentY = 2.0 * flatTopGaussianExponentY;
    xx = abs(Xrot-xorot);
    yy = abs(Yrot-yorot);
    F = gain * (...
        exp(-(xx/RcX).^exponentX) .* ...
        exp(-(yy/RcY).^exponentY));
end