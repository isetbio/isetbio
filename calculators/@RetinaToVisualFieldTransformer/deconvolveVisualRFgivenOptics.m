function retinalRFDoGparams = deconvolveVisualRFgivenOptics(obj, ...
    visualRFDoGparams, conesNumPooledByTheRFcenter, eccDegs, ...
    whichEye, testSubjectID, pupilDiameterMM, ...
    regularizationAlphaRange, varargin)
% Deconvolve a visual RF (DoG) assuming to retrieve the corresponding retinal RF
%
% Syntax:
%   retinalRFDoGparams = deconvolveVisualRFgivenOptics(obj, visualRFDoGparams, ...
%          conesNumPooledByTheRFcenter, eccDegs, ...
%          whichEye, testSubjectID, pupilDiameterMM, ...
%          regularizationAlphaRange, varargin)
%
% Description:
%    An object of the @cMosaic class contains parameters to create a
%    retinal cone mosaic and enable computation of cone excitations and
%
% Inputs:
%    None required.
%
% Outputs:
%    deconvolvedVisualRFdata
%
% Optional key/value pairs:
%    'wavelengthSupport'     - Vector. Wavalength support for computing the
%                                    Vlambda weighted optics
%
% History:
%    July 2025  NPC  Wrote it

    p = inputParser;
    p.addParameter('wavelengthSupport', 550 + 20*(-8:8), @isnumeric);
    p.parse(varargin{:});

    % Generate vLambda weighted PSF and OTF
    wavelengthSupport = 550 + 20*(-8:8);

    % Generate the cone mosaic at the given subject, eye and eccentricity
    cm = cMosaic(...
        'whichEye', whichEye, ...
        'sizeDegs', [0.5 0.5], ...
        'eccentricityDegs', eccDegs, ...
        'rodIntrusionAdjustedConeAperture', true, ...
        'customDegsToMMsConversionFunction', @RGCmodels.Watson.convert.rhoDegsToMMs, ...
        'customMMsToDegsConversionFunction', @RGCmodels.Watson.convert.rhoMMsToDegs, ...
        'wave', wavelengthSupport);

    % Generate Vlambda weighted PSF/OTF for the given subject, eye, pupil and eccentricity
    [theOTFData, thePSFData] = obj.vLambdaWeightedPSFandOTF(cm, testSubjectID, pupilDiameterMM);

    spatialSupportDegs = [thePSFData.supportX(:) thePSFData.supportY(:)]/60;

    % Generate a circularly-symmetric PSF
    circularSymmetryGenerationMode = 'average';
    circularSymmetryGenerationMode = 'bestResolution';
    %circularSymmetryGenerationMode = 'worseResolution';
    theCircularPSFData = thePSFData;
    theCircularPSFData.data = RetinaToVisualFieldTransformer.circularlySymmetricPSF(thePSFData.data, circularSymmetryGenerationMode); 
   
    if (~isempty(visualRFDoGparams.RcDegs))
        retinalRcDegs = [];
    else
        % Compute the visual Rc based on the retinal cone Rc, 
        % the number of cones in the RF center and their spacing and the PSF,
        % all for the current eccentricity

        % Sort cones according to their distance from the mosaic center
        coneDistancesFromMosaicCenter = sqrt(sum(bsxfun(@minus, cm.coneRFpositionsDegs, cm.eccentricityDegs).^2,2));
        [~,idx] = sort(coneDistancesFromMosaicCenter, 'ascend');

        % Estimate mean anatomical cone aperture in the mosaic'c center
        sourceConesIndices = idx(1:6);
        meanConeApertureDegsInMosaicCenter = mean(cm.coneApertureDiametersDegs(sourceConesIndices));
        anatomicalConeCharacteristicRadiusDegs = 0.204 * sqrt(2.0) * meanConeApertureDegsInMosaicCenter;

        % Compute the retinal RcDegs of the RF center composed of conesNumInRFcenter cones
        rfCenterPooledConeIndices = idx(1:conesNumPooledByTheRFcenter);
        conePosDegsRelativeToCenter = bsxfun(@minus, cm.coneRFpositionsDegs(rfCenterPooledConeIndices,:), cm.eccentricityDegs);
        [~, retinalRFcenterConeMap] = computeRetinalRFRcDegsFromItsPooledConeInputs(anatomicalConeCharacteristicRadiusDegs,...
            conePosDegsRelativeToCenter, spatialSupportDegs);

        % Convolve retinal cone map with the PSF
        visualRFcenterConeMap = conv2(theCircularPSFData.data, retinalRFcenterConeMap, 'same');

        % Fit a 2D Gaussian to the visually projected cone aperture and extract
        % the characteristic radius of that Gaussian
        [~,visualConeCharacteristicMinorMajorRadiiDegs] = RetinaToVisualFieldTransformer.fitGaussianEllipsoid(...
            thePSFData.supportX, thePSFData.supportY, visualRFcenterConeMap);

        visualRFDoGparams.RcDegs = min(visualConeCharacteristicMinorMajorRadiiDegs);
    end

    % Generate the visual RF (DoG with visualRFDoGparams
    visualDoGRFparamsVector = [...
        visualRFDoGparams.Kc ...
        visualRFDoGparams.RcDegs ...
        visualRFDoGparams.surroundToCenterRcRatio ...
        visualRFDoGparams.surroundToCenterIntegratedRatio ...
        ];
    RF2DData.visualRF = DoGRF(visualDoGRFparamsVector, spatialSupportDegs);
    RF2DData.visualRFcenterConeMap = visualRFcenterConeMap;
    RF2DData.retinalRFcenterConeMap = retinalRFcenterConeMap;

    % Deconvolve visual RF using the PSF data using the regularized filter
    % algorithm. The higher the alpha, the smoother the deconvolved image
    additiveNoisePower = 0;
    RF2DData.retinalRF = deconvreg(RF2DData.visualRF,theCircularPSFData.data, additiveNoisePower, regularizationAlphaRange);

    % Fit DoG model to the retinal RF to extract the corresponding params
    [retinalDoGRFparamsVector, RF2DData.fittedRetinalRF] = fitDoGModelToRF(...
        spatialSupportDegs, RF2DData.retinalRF, []);

    retinalRFDoGparams.Kc = retinalDoGRFparamsVector(1);
    retinalRFDoGparams.RcDegs = retinalDoGRFparamsVector(2);
    retinalRFDoGparams.surroundToCenterRcRatio = retinalDoGRFparamsVector(3);
    retinalRFDoGparams.surroundToCenterIntegratedRatio = retinalDoGRFparamsVector(4);

    % Convolved the retinal RF to see if we get back the visual RF
    RF2DData.visualRFcorrespondingToRetinalRF = conv2(RF2DData.retinalRF, theCircularPSFData.data, 'same');
    RF2DData.visualRFcorrespondingToFittedRetinalRF = conv2(RF2DData.fittedRetinalRF, theCircularPSFData.data, 'same');


    % Fit DoG model to the visual RF to extract the corresponding params
    [achievedVisualDoGRFparamsVector, RF2DData.fittedVisualRF] = fitDoGModelToRF(...
        spatialSupportDegs, RF2DData.visualRFcorrespondingToRetinalRF, visualRFDoGparams.RcDegs);

    achievedVisualRFDoGparams.Kc = achievedVisualDoGRFparamsVector(1);
    achievedVisualRFDoGparams.RcDegs = achievedVisualDoGRFparamsVector(2);
    achievedVisualRFDoGparams.surroundToCenterRcRatio = achievedVisualDoGRFparamsVector(3);
    achievedVisualRFDoGparams.surroundToCenterIntegratedRatio = achievedVisualDoGRFparamsVector(4);
    
    % Generate figure
    generateFigure(thePSFData, theCircularPSFData, RF2DData, ...
        visualRFDoGparams, retinalRFDoGparams, achievedVisualRFDoGparams, eccDegs, testSubjectID, 1);
end



function [theFittedRFDoGparams, theFittedRF] = fitDoGModelToRF(spatialSupportDegs, RF, RcDegs)
    
    %                             Kc   RcDegs   Rs/Rc integratedS/C
    DoGRFparams.initialValues = [  1    0.05    5     0.5];
    DoGRFparams.lowerBounds   = [ 1e-3  0.1/60  1     0.0];
    DoGRFparams.upperBounds   = [ 1e2   2       10     10];

    if (~isempty(RcDegs))
        DoGRFparams.initialValues(2) = RcDegs;
        DoGRFparams.lowerBounds(2) = RcDegs;
        DoGRFparams.upperBounds(2) = RcDegs;
    end

    % The optimization objective
    objective = @(p) sum((DoGRF(p, spatialSupportDegs) - RF).^2, 'all');

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
          'objective', objective, ...
          'x0', DoGRFparams.initialValues, ...
          'lb', DoGRFparams.lowerBounds, ...
          'ub', DoGRFparams.upperBounds, ...
          'options', options...
          );
      
     ms = MultiStart(...
          'Display', 'off', ...
          'StartPointsToRun','bounds-ineqs', ...  % run only initial points that are feasible with respect to bounds and inequality constraints.
          'UseParallel', true);
      
     % Run the multi-start
     multiStartsNum = 64;
     theFittedRFDoGparams = run(ms, problem, multiStartsNum);

     theFittedRF = DoGRF(theFittedRFDoGparams, spatialSupportDegs);
end


function generateFigure(thePSFData, theCircularPSFData, RF2DData, ...
    visualRFDoGparams, retinalRFDoGparams, achievedVisualRFDoGparams, eccDegs, testSubjectID, figNo)

    maxSpatialSupportDegs = 0.1;
    maxPSF = max([max(thePSFData.data(:)) max(theCircularPSFData.data(:))]);
    maxRF  = max([...
        max(abs(RF2DData.visualRF(:))) ...
        max(abs(RF2DData.retinalRF(:))) ...
        max(abs(RF2DData.visualRFcorrespondingToRetinalRF(:)))]);

    maxRFcenterConeMap = max([ ...
        max(RF2DData.visualRFcenterConeMap(:)) ...
        max(RF2DData.retinalRFcenterConeMap(:)) ...
        ]);

    hFig = figure(figNo); clf;
    set(hFig, 'Position', [15 15 1350 1000], 'Color', [1 1 1]);

    % The original PSF
    ax = subplot(3,4,1);
    plotPSF(ax, thePSFData, maxPSF, maxSpatialSupportDegs, eccDegs, testSubjectID);


    % The circular PSF
    ax = subplot(3,4,2);
    plotPSF(ax, theCircularPSFData, maxPSF, maxSpatialSupportDegs, eccDegs, testSubjectID);
    title('circularly-averaged PSF')

    % The retinal RF center cone map
    ax = subplot(3,4,3);
    plotRF(ax, thePSFData, RF2DData.retinalRFcenterConeMap, maxRFcenterConeMap, maxSpatialSupportDegs, 'retinal RF center');


    % The visual RF center cone map
    ax = subplot(3,4,4);
    plotRF(ax, thePSFData, RF2DData.visualRFcenterConeMap, maxRFcenterConeMap, maxSpatialSupportDegs, 'visual RF center');

    
    % The visual RF
    ax = subplot(3,4,5);
    plotRF(ax, thePSFData, RF2DData.visualRF, maxRF, maxSpatialSupportDegs, 'visual RF');


    % The retinal RF (via deconvolution of the visual RF
    ax = subplot(3,4,6);
    plotRF(ax, thePSFData, RF2DData.retinalRF, maxRF, maxSpatialSupportDegs, 'retinal RF (deconvolved visual RF)');


    % The fitted retinal RF 
    ax = subplot(3,4,7);
    plotRF(ax, thePSFData, RF2DData.fittedRetinalRF, maxRF, maxSpatialSupportDegs, 'fitted retinal RF (DoG model)');


    % The visual RF from the fitted retinal RF
    ax = subplot(3,4,8);
    plotRF(ax, thePSFData, RF2DData.visualRFcorrespondingToRetinalRF, maxRF, maxSpatialSupportDegs, 'visual RF (from retinal RF)');


    % Comparison of target visual and obtained visual RF
    ax = subplot(3,4,9);
    midPoint = (size(RF2DData.visualRF,1)-1)/2+1;
    plot(ax, thePSFData.supportX/60, RF2DData.visualRF(midPoint,:), 'k-', 'LineWidth', 1.5);
    hold(ax, 'on');
    plot(ax, thePSFData.supportX/60, RF2DData.visualRFcorrespondingToRetinalRF(midPoint,:), 'r--', 'LineWidth', 1.5);
    set(ax, 'XLim', maxSpatialSupportDegs*[-1 1], 'YLim', max([max(RF2DData.visualRF(:)) max(RF2DData.visualRFcorrespondingToRetinalRF(:))])*[-0.5 1], ...
        'XTick', -0.5:0.05:05, 'FontSize', 14);
    axis(ax, 'square');
    legend({'target', 'achieved'});
    title('target vs achieved visual RF')
    xlabel(ax,'degrees');

    
   
    % The DoG params
    ax = subplot(3,4,10);
    XYLims = [0.01 0.1];
    h1 = plot(ax,visualRFDoGparams.RcDegs, retinalRFDoGparams.RcDegs, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5]);
    hold(ax, 'on');
    h2 = plot(ax,achievedVisualRFDoGparams.RcDegs, retinalRFDoGparams.RcDegs, 'rs', 'MarkerSize', 8, 'MarkerFaceColor', [1 1 1]);
    
    plot(ax, XYLims, XYLims, 'k-');
    axis(ax, 'square');
    set(ax, 'XLim', XYLims, 'YLim', XYLims, ...
            'XTick', [0.001 0.003 0.01 0.03 0.1 0.3], 'YTick', [0.001 0.003 0.01 0.03 0.1 0.3], ...
            'YScale', 'log', 'XScale', 'log', 'FontSize', 14);
    grid(ax, 'on')
    legend(ax,[h1 h2], {'intended', 'achieved'});
    xlabel(ax, 'visual Rc (degs)');
    ylabel(ax, 'retinal Rc (degs)');

    ax = subplot(3,4,11);
    XYLims = [0.5 8];
    h1 = plot(ax,visualRFDoGparams.surroundToCenterRcRatio, retinalRFDoGparams.surroundToCenterRcRatio, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5]);
    hold(ax, 'on');
    h2 = plot(ax,achievedVisualRFDoGparams.surroundToCenterRcRatio, retinalRFDoGparams.surroundToCenterRcRatio, 'rs', 'MarkerSize', 8, 'MarkerFaceColor', [1 1 1]);
    legend(ax,[h1 h2], {'intended', 'achieved'});
    plot(ax, XYLims, XYLims, 'k-');
    axis(ax, 'square');
    set(ax, 'XLim', XYLims, 'YLim', XYLims, ...
            'XTick', 1:10, 'YTick', 1:10, 'FontSize', 14);
    grid(ax, 'on')
    xlabel(ax, 'visual Rs/Rc');
    ylabel(ax, 'retinal Rs/Rc');

    ax = subplot(3,4,12);
    XYLims = [0.1 10];
    h1 = plot(ax,visualRFDoGparams.surroundToCenterIntegratedRatio, retinalRFDoGparams.surroundToCenterIntegratedRatio, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5]);
    hold(ax, 'on');
    h2 = plot(ax,achievedVisualRFDoGparams.surroundToCenterIntegratedRatio, retinalRFDoGparams.surroundToCenterIntegratedRatio, 'rs', 'MarkerSize', 8, 'MarkerFaceColor', [1 1 1]);
    legend(ax,[h1 h2], {'intended', 'achieved'});
    plot(ax, XYLims, XYLims, 'k-');
    axis(ax, 'square');
    set(ax, 'XLim', XYLims, 'YLim', XYLims, ...
            'XTick', [0.1 0.3 1 3 10], 'YTick', [0.1 0.3 1 3 10], ...
            'YScale', 'log', 'XScale', 'log', 'FontSize', 14);
    grid(ax, 'on')
    xlabel(ax, 'visual int. S/C ratio');
    ylabel(ax, 'retinal int. S/C ratio');

end

function plotPSF(ax, thePSFData, maxPSF, maxSpatialSupportDegs, eccDegs, testSubjectID)
    psfZLevels = 0.05:0.1:0.95;
    contourf(ax,thePSFData.supportX/60, thePSFData.supportY/60, thePSFData.data/maxPSF, psfZLevels);
    hold on;
    midRow = (size(thePSFData.data,1)-1)/2+1;
    plot(ax, thePSFData.supportX/60, -maxSpatialSupportDegs*0.75 + 1.7*thePSFData.data(midRow,:)/maxPSF*maxSpatialSupportDegs, 'r-', 'LineWidth', 1.5);
    axis(ax,'image'); axis 'xy';
    set(ax, 'XLim', maxSpatialSupportDegs*[-1 1], 'YLim', maxSpatialSupportDegs*[-1 1], ...
        'XTick', -0.5:0.05:05, 'YTick', -0.5:0.05:05, 'CLim', [0 1], 'FontSize', 14);
    grid(ax, 'on');
    xlabel(ax,'degrees');
    title(ax, sprintf('PSF @(%2.0f,%2.0f) degs, subj: %d', eccDegs(1), eccDegs(2), testSubjectID));
    colormap(ax,brewermap(1024, 'greys'));
end


function plotRF(ax, thePSFData, RF, maxRF, maxSpatialSupportDegs, titleString)
    rfZLevels = -0.9:0.1:0.9;
    contourf(ax,thePSFData.supportX/60, thePSFData.supportY/60, RF/maxRF, rfZLevels);
    %imagesc(ax, thePSFData.supportX/60, thePSFData.supportY/60, RF/maxRF);
    hold on;
    midRow = (size(RF,1)-1)/2+1;
    plot(ax, thePSFData.supportX/60, -maxSpatialSupportDegs*0.5 + 1.4*RF(midRow,:)/maxRF*maxSpatialSupportDegs, 'r-', 'LineWidth', 1.5);
    plot(ax, thePSFData.supportX/60, -maxSpatialSupportDegs*0.5 + thePSFData.supportX*0, 'k-');
    axis(ax,'image'); axis 'xy';
    set(ax, 'XLim', maxSpatialSupportDegs*[-1 1], 'YLim', maxSpatialSupportDegs*[-1 1], ...
        'XTick', -0.5:0.05:05, 'YTick', -0.5:0.05:05, 'CLim', [-1 1], 'FontSize', 14);
    grid(ax, 'on');
    colormap(ax,brewermap(1024, '*RdBu'));
    xlabel(ax,'degrees');
    title(ax, titleString);
end

function RF2D = DoGRF(p, spatialSupportDegs) 
    [Xdegs,Ydegs] = meshgrid(spatialSupportDegs(:,1), spatialSupportDegs(:,2));
    Kc = p(1); RcDegs = p(2); RsDegs = RcDegs * p(3);
    Ks = Kc * p(4)/((RsDegs/RcDegs)^2);
    RF2D = Kc * exp(-(Xdegs/RcDegs).^2).*exp(-(Ydegs/RcDegs).^2) - Ks * exp(-(Xdegs/RsDegs).^2) .* exp(-(Ydegs/RsDegs).^2);

    rcToPixelResolutionRatio = RcDegs/(spatialSupportDegs(2,1)-spatialSupportDegs(1,1));
    if (rcToPixelResolutionRatio < 1)
       error('RcDegs is too small compared to PSF resolution');
    end 
end


function [RFcenterRcDegs, RF2D] = computeRetinalRFRcDegsFromItsPooledConeInputs(coneRcDegs, conePosDegs, spatialSupportDegs)
    
    [Xdegs,Ydegs] = meshgrid(spatialSupportDegs(:,1), spatialSupportDegs(:,2));

    conesNumPooledByRFcenter = size(conePosDegs,1);
    for iCone = 1:conesNumPooledByRFcenter
  
        theConeApertureRF = exp(-((Xdegs-conePosDegs(iCone,1))/coneRcDegs).^2).*exp(-((Ydegs-conePosDegs(iCone,2))/coneRcDegs).^2);
        
        if (iCone == 1)
            RF2D = theConeApertureRF;
        else
            RF2D = RF2D + theConeApertureRF;
        end
    end
    RF2D = RF2D / max(RF2D(:));

    centroidOfPooledCones = mean(conePosDegs,1);

    initialParams = [...
        max(RF2D(:)), ...
        centroidOfPooledCones(1), ...
        coneRcDegs * sqrt(conesNumPooledByRFcenter), ...
        centroidOfPooledCones(2), ...
        coneRcDegs * sqrt(conesNumPooledByRFcenter), ...
        0];

    lowerBounds = [...
        0, ...
        centroidOfPooledCones(1), ...
        coneRcDegs, ...
        centroidOfPooledCones(2), ...
        coneRcDegs, ...
        -180];

    upperBounds = [...
        10, ...
        centroidOfPooledCones(1), ...
        coneRcDegs*conesNumPooledByRFcenter, ...
        centroidOfPooledCones(2), ...
        coneRcDegs*conesNumPooledByRFcenter, ...
        360];

    [theFittedGaussianCharacteristicRadiusDegs, theFittedGaussianEllpsoid, XYcenter, XRange, YRange, theFittedGaussianCharacteristicRadiiDegs] = ...
        fitGaussianToPooledConeApertures(spatialSupportDegs(:,1), spatialSupportDegs(:,2), RF2D, ...
        initialParams, lowerBounds, upperBounds);

    % Return the smallest of the 2 characteristic radii
    RFcenterRcDegs = max(theFittedGaussianCharacteristicRadiiDegs); %theFittedGaussianCharacteristicRadiusDegs; % min(theFittedGaussianCharacteristicRadiiDegs);

    debugFitGaussianToPooledConeApertures = false;
    if (debugFitGaussianToPooledConeApertures)
        figure(234); clf;
        zLevels = 0.05:0.1:1;
        subplot(1,2,1);
        contourf(spatialSupportDegs(:,1), spatialSupportDegs(:,2), RF2D, zLevels);
        hold on;
        contour(spatialSupportDegs(:,1), spatialSupportDegs(:,2), theFittedGaussianEllpsoid, zLevels, 'LineColor', 'r');
        axis 'image'; axis 'xy'
        set(gca, 'XLim', 0.05*[-1 1], 'YLim', 0.05*[-1 1], 'CLim', [0 1])
    
        colormap(gray(1024))
        pause
    end

end

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
