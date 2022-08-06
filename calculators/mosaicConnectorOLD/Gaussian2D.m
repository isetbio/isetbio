function Gaussian2D
% Set up mesh
xmin = -100;
xmax = 100;
ymin = -100;
ymax = 100;
spacing = 0.1;
xvals = xmin:spacing:xmax+spacing/2;
yvals = ymin:spacing:ymax+spacing/2;
[x,y] = meshgrid(xvals, yvals);
xydata(:,:,1) = x;
xydata(:,:,2) = y;

% parameters for the gaussian
params(1) = 1;  % amplitude
params(2) = 20; % x0
params(3) = 40; % y0
params(4) = 20; % xSigma
params(5) = 2; % ySigma
params(6) = 40; % rotation (degs)
f = elliptical2DGaussian(params,xydata);
idx = find(f >= exp(-1));
f = f / sum(f(:));

f1e = f(idx);
[sum(f1e(:)) sum(f(:))]
pause

% Compute the filter and display it
figure(1)
contour(x, y, f);

[fittedParams,resnorm,residual,exitflag] = lsqcurvefit(@elliptical2DGaussian,initialParams,xydata,theRF,lb,ub);


 
end


function f = elliptical2DGaussian(params,xydata)
    amplitude = params(1);
    xCenter = params(2);
    yCenter = params(3);
    xSigma = params(4);
    ySigma = params(5);
    theta = params(6);

    xx = xydata(:,:,1)*cosd(theta) - xydata(:,:,2)*sind(theta);
    yy = xydata(:,:,1)*sind(theta) + xydata(:,:,2)*cosd(theta);
    xx0 = xCenter*cosd(theta) - yCenter*sind(theta);
    yy0 = xCenter*sind(theta) + yCenter*cosd(theta);

    f = amplitude * exp( -((xx-xx0).^2/(xSigma^2) + (yy-yy0).^2/(ySigma^2)) );
end


