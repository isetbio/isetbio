% Method to fit a 2D Gaussian to the visually projected cone aperture
function [theFittedGaussianCharacteristicRadiusDegs, theFittedGaussianEllpsoid, XYcenter, XRange, YRange, theFittedGaussianCharacteristicRadiiDegs] = ...
    fitGaussianToPooledConeApertures(supportX, supportY, theRF, initialParams, lowerBounds, upperBounds)

    [X,Y] = meshgrid(supportX, supportY);
    xydata(:,:,1) = X;
    xydata(:,:,2) = Y;


    theRF = theRF / max(theRF(:));

    % Form parameter vector: [gain, xo, RcX, yo, RcY, rotationAngleDegs]
    p0 = initialParams;

    % Form lower and upper value vectors
    lb = lowerBounds;
    ub = upperBounds;

    % Do the fitting
    [fittedParams,resnorm,residual,exitflag] = lsqcurvefit(@RetinaToVisualFieldTransformer.gaussian2D,p0,xydata,theRF,lb,ub);

    xo = fittedParams(2);
    yo = fittedParams(4);
    RcX = fittedParams(3);
    RcY = fittedParams(5);
    XRange = max([RcX RcY])*[-3 3];
    YRange = max([RcX RcY])*[-3 3];
    XYcenter = [xo yo]; 

    % Compute the fitted 2D Gaussian
    theFittedGaussianEllpsoid = RetinaToVisualFieldTransformer.gaussian2D(fittedParams,xydata);

    % Compute the fitted Gaussian Rc in degs
    theFittedGaussianCharacteristicRadiusDegs = (sqrt(RcX^2+RcY^2)/sqrt(2));

    theFittedGaussianCharacteristicRadiiDegs = [RcX RcY];
end