function testLineVsPointSpread
    % Spatial support
    N = 1024;
    x = linspace(-5,5,N);deltaX = x(2)-x(1);
    [X,Y] = meshgrid(x,x);
    
    % Center characteristic radius and point sensitivity
    rc = 0.2; kc = 1;
    
    % Surround characteristic radius and point sensitivity
    rs = rc*5; ks = kc/50;
    
    % Classic DoG model point weighting function
    centerPointWeightingFunction2D = kc * exp(-( (X/rc).^2 + (Y/rc).^2 ));
    surroundPointWeightingFunction2D = ks * exp(-( (X/rs).^2 + (Y/rs).^2 ));
    pointWeightingFunction2D = centerPointWeightingFunction2D  - surroundPointWeightingFunction2D;
    
    % Analytically-derived line weighting function
    lineWeightingFunctionCenter = kc * rc * sqrt(pi)*exp(-(x/rc).^2);
    lineWeightingFunctionSurround = ks*rs*sqrt(pi)*exp(-(x/rs).^2);
    lineWeightingFunction = lineWeightingFunctionCenter - lineWeightingFunctionSurround;
    
    % Computationally-derived line weighting function
    lineWeightingFunctionFromIntegration = squeeze(sum(pointWeightingFunction2D, 2))*deltaX;
    
    % Contrast sensitivity function (1D gratings)
    sf = logspace(log10(0.01), log10(5), 24);
    K = 1;
    CSF = K * pi * (kc * rc^2*exp(-(pi*rc*sf).^2) - ks * rs^2*exp(-(pi*rs*sf).^2));
    CSF = CSF/max(CSF);
    
    % Computationally-derived CSF
    FTmag = fftshift(abs(fft(lineWeightingFunction)));
    FTmag = FTmag / max(FTmag);
    sfMax = 1/(2*deltaX); deltaSF = sfMax/(N/2);
    FTsupport = linspace(-sfMax, sfMax-deltaSF, N);

    figure(1); clf;
    subplot(2,2,1);
    % Extract a slice though the point weighting function
    midPoint = floor(size(pointWeightingFunction2D,1)/2);
    pointWeightingFunctionSlice = squeeze(pointWeightingFunction2D(midPoint,:));
    pointWeightingFunctionSlice = pointWeightingFunctionSlice / max(pointWeightingFunctionSlice);
    
    centerPointWeightingFunctionSlice = squeeze(centerPointWeightingFunction2D(midPoint,:));
    surroundPointWeightingFunctionSlice = squeeze(surroundPointWeightingFunction2D(midPoint,:));
    
    plot(x, pointWeightingFunctionSlice, 'ks-'); hold on;
    plot(x, lineWeightingFunctionFromIntegration, 'ro-');
    plot(x, lineWeightingFunction, 'b.');
    set(gca, 'XLim', 2*[-1 1]);
    legend({'point weighting', 'line weighting (int)', 'lineWeightingFunction'});
    
    subplot(2,2,2);
    plot(x, lineWeightingFunction/max(lineWeightingFunction), 'b.');
    set(gca, 'XLim', 2*[-1 1]);
    
    subplot(2,2,3)
    plot(x, centerPointWeightingFunctionSlice, 'r-'); hold on;
    plot(x, surroundPointWeightingFunctionSlice, 'b-');
    plot(x, lineWeightingFunctionCenter, 'rs-');
    plot(x, lineWeightingFunctionSurround, 'bs-');
    
    subplot(2,2,4);
    plot(sf, CSF, 'ks-'); hold on;
    plot(FTsupport, FTmag, 'ro');
    set(gca, 'XScale' ,'log', 'YScale', 'linear', 'XLim', [0.01 5]);
    
end

