% Method to fit a 2D Gaussian to the visually projected cone aperture
function [theFittedGaussianCharacteristicRadiusDegs, ...
    theFittedGaussianCharacteristicRadiiDegs, ...
    theFittedGaussianRotationDegs, ...
    theFittedGaussianFlatTopExponent, ...
    theFittedGaussianEllpsoid, XYcenter, XRange, YRange] = ...
    fitGaussianEllipsoid(supportX, supportY, theRF, varargin)

    p = inputParser;
    p.addParameter('flatTopGaussian', false, @islogical);
    p.addParameter('forcedOrientationDegs', [], @(x)(isempty(x) || isscalar(x)));
    p.parse(varargin{:});
    flatTopGaussian = p.Results.flatTopGaussian;
    forcedOrientationDegs = p.Results.forcedOrientationDegs;

    [X,Y] = meshgrid(supportX, supportY);
    xydata(:,:,1) = X;
    xydata(:,:,2) = Y;

    maxRF = max(theRF(:));
    theRF = theRF / maxRF;
    [theCentroid, RcX, RcY, theRotationAngle] = RetinaToVisualFieldTransformer.estimateGeometry(supportX, supportY, theRF);

    if (~isempty(forcedOrientationDegs))
        theRotationAngle = forcedOrientationDegs;
    end

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

    if (~isempty(forcedOrientationDegs))
        lb(end) = forcedOrientationDegs;
        ub(end) = forcedOrientationDegs;
    end

    if (flatTopGaussian)
        p0(numel(p0)+1) = 0.6;
        lb(numel(lb)+1) = 0.25;
        ub(numel(ub)+1) = 1.0;
    else
        p0(numel(p0)+1) = 1.0;
        lb(numel(lb)+1) = 1.0;
        ub(numel(ub)+1) = 1.0;
    end

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
    theFittedGaussianRotationDegs = fittedParams(6);
    theFittedGaussianFlatTopExponent = fittedParams(7);
    
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
    flatTopGaussianExponent = params(7);

    % Apply axes rotation
    Xrot = X * cosd(rotationAngle) -  Y*sind(rotationAngle);
    Yrot = X * sind(rotationAngle) +  Y*cosd(rotationAngle);
    xorot = xo * cosd(rotationAngle) -  yo*sind(rotationAngle);
    yorot = xo * sind(rotationAngle) +  yo*cosd(rotationAngle);

    % Compute 2D Gaussian
    F = gain * (exp(-((Xrot-xorot)/RcX).^2) .* exp(-((Yrot-yorot)/RcY).^2)).^flatTopGaussianExponent;
end