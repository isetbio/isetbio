function dStruct = estimateConeCharacteristicRadiusInVisualSpace(obj,...
            whichEye, eccDegs, testSubjectID, pupilDiameterMM,  dataFileName, varargin)
% Estimate the characteristic radius of the mean cone at the center of a
% cone mosaic centered at a given (x,y) eccentricity in degrees
% as projected on to visual space using optics at the same eccentricity

    % Parse input
    p = inputParser;
    p.addParameter('anatomicalConeCharacteristicRadiusDegs', [], @(x)(isempty(x)||isscalar(x)));
    p.addParameter('hFig', []);
    p.addParameter('videoOBJ', []);
    p.parse(varargin{:});

    anatomicalConeCharacteristicRadiusDegs = p.Results.anatomicalConeCharacteristicRadiusDegs;
    hFig = p.Results.hFig;
    videoOBJ = p.Results.videoOBJ;

    % Generate the cone mosaic at the given eye and eccentricity
    cm = cMosaic(...
        'whichEye', whichEye, ...
        'sizeDegs', [1 1], ...
        'eccentricityDegs', eccDegs, ...
        'customDegsToMMsConversionFunction', @RGCmodels.Watson.convert.rhoDegsToMMs, ...
        'customMMsToDegsConversionFunction', @RGCmodels.Watson.convert.rhoMMsToDegs, ...
        'wave', 550);

    % How many cones we have in this [1x1] deg patch
    conesNumInRetinalPatch = numel(cm.coneTypes);

    % If an anatomicalConeCharacteristicRadiusDegs) is not specified,
    % compute it here from the center of the mosaic
    if (isempty(anatomicalConeCharacteristicRadiusDegs))
        % If we have less than 6 cones, we will not do any analysis
        if (conesNumInRetinalPatch < 6)
            fprintf(2, 'Less than 6 cones at this eccentricity, skipping computation of cone aperture in visual space.\n')
            
            % Return struct
            dStruct.conesNumInRetinalPatch = conesNumInRetinalPatch;
            dStruct.anatomicalConeCharacteristicRadiusDegs = nan;
            dStruct.visualConeCharacteristicRadiusDegs = nan;
            return;
        end

        % Estimate mean anatomical cone aperture in the mosaic'c center
        coneDistancesFromMosaicCenter = sqrt(sum(bsxfun(@minus, cm.coneRFpositionsDegs, cm.eccentricityDegs).^2,2));
        [~,idx] = sort(coneDistancesFromMosaicCenter);
        meanConeSpacingDegsInMosaicCenter = mean(cm.coneRFspacingsDegs(idx(1:6)));
        anatomicalConeCharacteristicRadiusDegs = 0.204 * sqrt(2.0) * meanConeSpacingDegsInMosaicCenter;
    end


    % Generate optics for this eye, eccentricity, subject, and pupil size
    switch (obj.ZernikeDataBase)
        % Artal
        case RetinaToVisualFieldTransformer.Artal
            subtractCentralRefraction = ArtalOptics.constants.subjectRequiresCentralRefractionCorrection(...
                cm.whichEye, testSubjectID);
        % Polans
        case RetinaToVisualFieldTransformer.Polans
            subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(...
                cm.whichEye, testSubjectID);
    end

    [oiEnsemble, psfEnsemble] = cm.oiEnsembleGenerate(cm.eccentricityDegs, ...
                    'zernikeDataBase', obj.ZernikeDataBase, ...
                    'subjectID', testSubjectID, ...
                    'pupilDiameterMM', pupilDiameterMM, ...
                    'zeroCenterPSF', true, ...
                    'subtractCentralRefraction', subtractCentralRefraction, ...
                    'wavefrontSpatialSamples', 701, ...
                    'warningInsteadOfErrorForBadZernikeCoeffs', true);

    if (isempty(oiEnsemble))
        fprintf(2, 'Could not generate optics at this eccentricity.\n');
        % Return struct
        dStruct.conesNumInRetinalPatch = conesNumInRetinalPatch;
        dStruct.anatomicalConeCharacteristicRadiusDegs = anatomicalConeCharacteristicRadiusDegs;
        dStruct.visualConeCharacteristicRadiusDegs = nan;
        return;
    end

    % Extract the PSF
    thePSFData = psfEnsemble{1};
    targetWavelength = cm.wave(1);
    [~, wIndex] = min(abs(thePSFData.supportWavelength-targetWavelength));
    thePSFData.data = squeeze(thePSFData.data(:,:,wIndex));


    p = getpref('ISETMacaque');
    pdfFileName = fullfile(p.generatedDataDir, 'coneApertureBackingOut', ...
            sprintf('%s_ecc_%2.0f_%2.0f.pdf',dataFileName, cm.eccentricityDegs(1), cm.eccentricityDegs(2)));

    % Compute the visually-projected cone aperture given this PSF
    visualConeCharacteristicRadiusDegs = analyzeEffectOfPSFonConeAperture(...
                    anatomicalConeCharacteristicRadiusDegs, thePSFData, ...
                    hFig, videoOBJ, pdfFileName);

    % Return struct
    dStruct.conesNumInRetinalPatch = conesNumInRetinalPatch;
    dStruct.anatomicalConeCharacteristicRadiusDegs = anatomicalConeCharacteristicRadiusDegs;
    dStruct.visualConeCharacteristicRadiusDegs = visualConeCharacteristicRadiusDegs;
end


function visualConeCharacteristicRadiusDegs = analyzeEffectOfPSFonConeAperture(...
     anatomicalConeCharacteristicRadiusDegs, thePSFData, ...
     hFig, videoOBJ, pdfFileName)

    [Xarcmin, Yarcmin] = meshgrid(thePSFData.supportX, thePSFData.supportY);

    [theCentroid, RcX, RcY, theRotationAngle] = ...
        estimatePSFgeometry(thePSFData.supportX, thePSFData.supportY,thePSFData.data);

    % Generate anatomical cone aperture (a Gaussian with Rc)
    anatomicalConeAperture = exp(-(((Xarcmin-theCentroid(1))/60)/anatomicalConeCharacteristicRadiusDegs).^2) .* ...
                             exp(-(((Yarcmin-theCentroid(2))/60)/anatomicalConeCharacteristicRadiusDegs).^2);

    % Convolve cone aperture with the PSF
    theVisuallyProjectedConeAperture = conv2(thePSFData.data, anatomicalConeAperture, 'same');
    theVisuallyProjectedConeAperture = theVisuallyProjectedConeAperture  / max(theVisuallyProjectedConeAperture(:));

    % Fit a 2D Gaussian to the visually projected cone aperture and extract
    % the characteristic radius of that Gaussian
    [visualConeCharacteristicRadiusDegs, ...
     theVisuallyProjectedConeApertureFittedGaussian, XYcenter, XRange, YRange] = ...
        fit2DGaussianToVisuallyProjectedConeAperture(thePSFData.supportX, thePSFData.supportY, theVisuallyProjectedConeAperture);

    if (isempty(hFig))   
        return;
    end

    XLims = XYcenter(1) + 7*[-1 1];  % XLims
    YLims = XYcenter(2) + 7*[-1 1];  % YLims
    plotAnalysis(anatomicalConeAperture, thePSFData,...
        theVisuallyProjectedConeAperture,...
        theVisuallyProjectedConeApertureFittedGaussian, ...
        XLims, YLims, videoOBJ, pdfFileName)
end

function plotAnalysis(anatomicalConeAperture, thePSFData, ...
    theVisuallyProjectedConeAperture, ...
    theVisuallyProjectedConeApertureFittedGaussian, ...
    XLims, YLims, videoOBJ, pdfFileName)
    
    xRange = XLims(2)-XLims(1);
    yRange = YLims(2)-YLims(1);
    xyRange = max([xRange yRange]);
    if (xyRange < 10)
        xTick = -100:1:100;
    elseif (xyRange < 30)
        xTick = -100:2:100;
    elseif (xyRange < 50)
        xTick = -100:5:100;
    elseif (xyRange < 100)
        xTick = -200:10:200;
    else
        xTick = -400:20:400;
    end

    hFig = figure(10); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1600 400]);
    
    
    % Plot here
    cLUT= brewermap(1024, 'blues');
    zLevels = 0.05:0.15:1.0;

    % The cone aperture
    ax = subplot(1,4,1);
    contourf(ax,thePSFData.supportX, thePSFData.supportY, anatomicalConeAperture, zLevels);
    set(ax, 'XLim', XLims, 'YLim', YLims, 'FontSize', 14, 'Color', squeeze(cLUT(1,:)), 'XTick', xTick, 'YTick', xTick);
    axis(ax,'xy'); axis(ax, 'square'); 
    grid(ax, 'on');
    xlabel(ax,'arc min');
    xtickangle(ax, 0);
    title(ax,'cone aperture');
    
    ax = subplot(1,4,2);
    contourf(ax,thePSFData.supportX, thePSFData.supportY, thePSFData.data/max(thePSFData.data(:)), zLevels);
    set(ax, 'XLim', XLims, 'YLim', YLims, 'FontSize', 14, 'Color', squeeze(cLUT(1,:)), 'XTick', xTick, 'YTick', xTick);
    axis(ax,'xy'); axis(ax, 'square'); 
    grid(ax, 'on');
    xlabel(ax,'arc min');
    xtickangle(ax, 0);
    title(ax,'point spread function');

    ax = subplot(1,4,3);
    contourf(ax,thePSFData.supportX, thePSFData.supportY, theVisuallyProjectedConeAperture, zLevels);
    set(ax, 'XLim', XLims, 'YLim', YLims, 'FontSize', 14, 'Color', squeeze(cLUT(1,:)), 'XTick', xTick, 'YTick', xTick);
    axis(ax,'xy'); axis(ax, 'square'); 
    grid(ax, 'on');
    xlabel(ax,'arc min');
    xtickangle(ax, 0);
    title(ax,sprintf('visually projected cone aperture\n conv(coneAperture, PSF)'));
    
    ax = subplot(1,4,4);
    contourf(ax,thePSFData.supportX, thePSFData.supportY, theVisuallyProjectedConeApertureFittedGaussian, zLevels);
    set(ax, 'XLim', XLims, 'YLim', YLims, 'FontSize', 14, 'Color', squeeze(cLUT(1,:)), 'XTick', xTick, 'YTick', xTick);
    axis(ax,'xy'); axis(ax, 'square'); 
    grid(ax, 'on');
    xlabel(ax,'arc min');
    xtickangle(ax, 0);
    title(ax,'fitted Gaussian ellipsoid');
    
    colormap(cLUT);

    drawnow;
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);

    if (~isempty(videoOBJ))
        videoOBJ.writeVideo(getframe(hFig));
    end

end


function [theCentroid, RcX, RcY, theRotationAngle] = estimatePSFgeometry(supportX, supportY, theVisualConeAperture)
    % Compute orientation, centroid, and major/minor axis lengths
    binaryImage = theVisualConeAperture;
    m1 = min(binaryImage(:));
    m2 = max(binaryImage(:));
    binaryImage = imbinarize((binaryImage - m1)/(m2-m1));
    s = regionprops('table', binaryImage,'Orientation', 'Centroid', 'MinorAxisLength', 'MajorAxisLength');
    theCentroid = s.Centroid(1,:);
    theMinorAxisLength = s.MinorAxisLength(1);
    theMajorAxisLength = s.MajorAxisLength(1);
    theRotationAngle = s.Orientation(1);

    % The computed centroid and axis lengths are in pixels. Convert them to units of spatial support
    % to serve as initial Gaussian parameter values
    xx = [round(theCentroid(1)-theMajorAxisLength*0.5) round(theCentroid(1)+theMajorAxisLength*0.5)];
    yy = [round(theCentroid(2)-theMinorAxisLength*0.5) round(theCentroid(2)+theMinorAxisLength*0.5)];

    xx(1) = max([1 xx(1)]);
    xx(2) = min([numel(supportX) xx(2)]);

    yy(1) = max([1 yy(1)]);
    yy(2) = min([numel(supportY) yy(2)]);

    RcY = (supportY(yy(2)) - supportY(yy(1)))/5.0;
    RcX = (supportX(xx(2)) - supportX(xx(1)))/5.0;
    theCentroid(1) = supportX(round(theCentroid(1)));
    theCentroid(2) = supportY(round(theCentroid(2)));
end


% Method to fit a 2D Gaussian to the visually projected cone aperture
function [theFittedGaussianCharacteristicRadiusDegs, theFitted2DGaussian, XYcenter, XRange, YRange] = ...
    fit2DGaussianToVisuallyProjectedConeAperture(supportX, supportY, theVisualConeAperture)

    [X,Y] = meshgrid(supportX, supportY);
    xydata(:,:,1) = X;
    xydata(:,:,2) = Y;

    theVisualConeAperture = theVisualConeAperture / max(theVisualConeAperture(:));
    [theCentroid, RcX, RcY, theRotationAngle] = estimatePSFgeometry(supportX, supportY, theVisualConeAperture);

    % Form parameter vector: [gain, xo, RcX, yo, RcY, rotationAngleDegs]
    p0 = [...
        max(theVisualConeAperture(:)), ...
        theCentroid(1), ...
        RcX, ...
        theCentroid(2), ...
        RcY, ...
        theRotationAngle];

    % Form lower and upper value vectors
    lb = [ 0 min(supportX) 0*(max(supportX)-min(supportX))    min(supportY) 0*(max(supportY)-min(supportY))  theRotationAngle-90];
    ub = [ 1 max(supportX)    max(supportX)-min(supportX)     max(supportY)    max(supportY)-min(supportY)   theRotationAngle+90];

    % Do the fitting
    [fittedParams,resnorm,residual,exitflag] = lsqcurvefit(@gaussian2D,p0,xydata,theVisualConeAperture,lb,ub);

    xo = fittedParams(2);
    yo = fittedParams(4);
    RcX = fittedParams(3);
    RcY = fittedParams(5);
    XRange = max([RcX RcY])*[-3 3];
    YRange = max([RcX RcY])*[-3 3];
    XYcenter = [xo yo]; 

    % Compute the fitted 2D Gaussian
    theFitted2DGaussian = gaussian2D(fittedParams,xydata);

    % Compute the fitted Gaussian Rc in degs
    theFittedGaussianCharacteristicRadiusDegs = sqrt(RcX^2+RcY^2)/60;
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