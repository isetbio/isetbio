% Method to fit a 2D Gaussian to the visually projected cone aperture
function [theFittedGaussianCharacteristicRadiusDegs, ...
    theFittedGaussianCharacteristicRadiiDegs, ...
    theFittedGaussianEllpsoid, XYcenter, XRange, YRange] = ...
    fitGaussianEllipsoid(supportX, supportY, theRF)

    [X,Y] = meshgrid(supportX, supportY);
    xydata(:,:,1) = X;
    xydata(:,:,2) = Y;

    maxRF = max(theRF(:));
    theRF = theRF / maxRF;
    [theCentroid, RcX, RcY, theRotationAngle] = RetinaToVisualFieldTransformer.estimateGeometry(supportX, supportY, theRF);

    % Form parameter vector: [gain, xo, RcX, yo, RcY, rotationAngleDegs]
    p0 = [...
        max(theRF(:)), ...
        theCentroid(1), ...
        RcX, ...
        theCentroid(2), ...
        RcY, ...
        theRotationAngle];

    % Form lower and upper value vectors
    lb = [ 0 min(supportX) 0*(max(supportX)-min(supportX))    min(supportY) 0*(max(supportY)-min(supportY))  theRotationAngle-90];
    ub = [ 1 max(supportX)    max(supportX)-min(supportX)     max(supportY)    max(supportY)-min(supportY)   theRotationAngle+90];

    % Do the fitting
    [fittedParams,resnorm,residual,exitflag] = lsqcurvefit(@gaussian2D,p0,xydata,theRF,lb,ub);

    xo = fittedParams(2);
    yo = fittedParams(4);
    RcX = fittedParams(3);
    RcY = fittedParams(5);
    XRange = max([RcX RcY])*[-3 3];
    YRange = max([RcX RcY])*[-3 3];
    XYcenter = [xo yo]; 

    % Compute the fitted 2D Gaussian
    theFittedGaussianEllpsoid = gaussian2D(fittedParams,xydata);

    % Compute the fitted Gaussian Rc in degs
    theFittedGaussianCharacteristicRadiusDegs = (sqrt(RcX^2+RcY^2)/sqrt(2));

    theFittedGaussianCharacteristicRadiiDegs = [RcX RcY];
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
    RcY = params(5);
    rotationAngle = params(6);

    % Apply axes rotation
    Xrot = X * cosd(rotationAngle) -  Y*sind(rotationAngle);
    Yrot = X * sind(rotationAngle) +  Y*cosd(rotationAngle);
    xorot = xo * cosd(rotationAngle) -  yo*sind(rotationAngle);
    yorot = xo * sind(rotationAngle) +  yo*cosd(rotationAngle);

    % Compute 2D Gaussian
    F = gain * exp(-((Xrot-xorot)/RcX).^2) .* exp(-((Yrot-yorot)/RcY).^2);
end