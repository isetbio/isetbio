function [ellipsoidParams, RFmapHR, xSupportHR, ySupportHR] = ...
    GaussianEllipsoidTo2DRFmap(xSupport, ySupport, theRFmap, nSamplesHighRes, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('fixedCentroid', [], @(x)(isempty(x)||(numel(x)==2)));
    p.parse(varargin{:});
    fixedCentroid = p.Results.fixedCentroid;

    options = optimoptions(@lsqcurvefit,'Algorithm','levenberg-marquardt','MaxIterations',1500, 'UseParallel', true);

	[X,Y] = meshgrid(xSupport, ySupport);
	spatialSupport = zeros(size(X,1),size(Y,2),2);
	spatialSupport(:,:,1) = X;
	spatialSupport(:,:,2) = Y;
	xyRange = (xSupport(end)-xSupport(1));

	maxSigmaX = (max(xSupport)/2)^2;
	maxSigmaY = (max(ySupport)/2)^2;
	%               [gain, xo,            xSigma,                 yo,           ySigma,            orientation]
	initialParams = [1     0              0.1*maxSigmaX            0             0.1*maxSigmaY     0.0];
    lowerBounds   = [0     min(xSupport)  0                       min(ySupport)  0                 -pi/2];
    upperBounds   = [2     max(xSupport)  maxSigmaX               max(ySupport)  maxSigmaY   		pi/2];

    rangeX = max(xSupport) - min(xSupport);
    rangeY = max(ySupport) - min(ySupport);

    if (~isempty(fixedCentroid))
    	% xo
        initialParams(2) = fixedCentroid(1);
        lowerBounds(2) = initialParams(2)-0.1*rangeX;
        upperBounds(2) = initialParams(2)+0.1*rangeX;
        % yo
        initialParams(4) = fixedCentroid(2);
        lowerBounds(4) = initialParams(4)-0.1*rangeY;
        upperBounds(4) = initialParams(4)+0.1*rangeY;
    end

    [bestFitParams,resnorm,residual,exitflag] = ...
    	lsqcurvefit(@GaussianEllipsoid, initialParams, spatialSupport, theRFmap, ...
        lowerBounds, upperBounds, options);

    if (exitflag == 1) && (~isempty(fixedCentroid))
    	% Fit failed. Try again with truly fixedCentroid (zero tolerance)
    	fprintf('----> Fit failed. Trying fitting once more, now with zero tolerance in the fixedCentroid.\n')
    	% xo
        initialParams(2) = fixedCentroid(1);
        lowerBounds(2) = initialParams(2);
        upperBounds(2) = initialParams(2);
        % yo
        initialParams(4) = fixedCentroid(2);
        lowerBounds(4) = initialParams(4);
        upperBounds(4) = initialParams(4);

        [bestFitParams,resnorm,residual,exitflag] = ...
    		lsqcurvefit(@GaussianEllipsoid, initialParams, spatialSupport, theRFmap, ...
        	lowerBounds, upperBounds, options);
    end

    if (exitflag == 1)
    	fprintf('----> Fit failed !!!\n');
    end


    ellipsoidParams.gain = bestFitParams(1);
	ellipsoidParams.x0 = bestFitParams(2);
	ellipsoidParams.xSigma = bestFitParams(3);
	ellipsoidParams.y0 = bestFitParams(4);
	ellipsoidParams.ySigma = bestFitParams(5);
	ellipsoidParams.rotationDegs = bestFitParams(6)/pi*180;


	if (~isempty(nSamplesHighRes))
		% Generate high resolution Gaussian Ellipsoid
	    nSamplesHR = 512;
	    xSupportHR = linspace(xSupport(1), xSupport(end), nSamplesHR);
	    ySupportHR = linspace(ySupport(1), ySupport(end), nSamplesHR);

	    [XHR,YHR] = meshgrid(xSupportHR, ySupportHR); 
	    xdataHR = zeros(nSamplesHR,nSamplesHR,2);
		spatialSupportHR(:,:,1) = XHR;
		spatialSupportHR(:,:,2) = YHR;

		RFmapHR = GaussianEllipsoid(bestFitParams, spatialSupportHR);
	else
		RFmapHR = [];
		xSupportHR = []; 
		ySupportHR = [];
	end
end

function F = GaussianEllipsoid(params, xdata)	
	gain = params(1);
	x0 = params(2);    
 	y0 = params(4);
	xSigma = params(3); 
	ySigma = params(5);
	theRotation = params(6);

	cosPhi = cos(theRotation);
	sinPhi = sin(theRotation);

	theRotationMatrix = [...
		cosPhi  -sinPhi; ...
		sinPhi   cosPhi];

	X = xdata(:,:,1); Y = xdata(:,:,2);
	n = size(X,1); m = size(X,2);

	XYrot = [X(:) Y(:)] * theRotationMatrix;
	xy0rot = [x0 y0] * theRotationMatrix;

	X = XYrot(:,1) - xy0rot(1);
	Y = XYrot(:,2) - xy0rot(2);

	F = gain * exp(-0.5*(X/xSigma).^2) .* exp(-0.5*(Y/ySigma).^2);
	F = reshape(F, [n m]);
end

