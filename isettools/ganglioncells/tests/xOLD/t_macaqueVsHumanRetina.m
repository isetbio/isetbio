function t_macaqueVsHumanRetina

    % Function handle to derive linear distance from angular distance on the macaque retina (222 microns/deg)
    micronsPerDegreeInMacaqueRetina = 222;
    macaqueAngularDistanceDegstoLinearDistanceMMs = @(angularDegs)(angularDegs * micronsPerDegreeInMacaqueRetina * 1e-3);

    % Load macaque cone density data from packer
    dMacaque = RGCmodels.PackerDacey.coneSpacingAndDensityAlongHorizontalMeridian();

    % Retrieve the angular eccentricity (temporal meridian)
    eccDegsOnMacaqueRetina = dMacaque.temporalMeridianData.eccDegs;
    eccMMOnMacaqueRetina = dMacaque.temporalMeridianData.eccMM;

    % Assume human and macaque retinas are equivalent at identical linear
    % ecccentricities
    eccMMOnHumanRetina = eccMMOnMacaqueRetina;

    % Derive angular eccentricity in human retina
    eccDegsOnHumanRetina = RGCmodels.Watson.convert.rhoMMsToDegs(eccMMOnHumanRetina);

    % Obtain human cone density data along the temporal meridian
    useParfor = true;
    WatsonMeridianName = RGCmodels.Watson.convert.retinalMeridianNameToVisualMeridianName('temporal meridian', 'right eye');
    [coneSpacingDegs, ~, ~, coneDensityMMs2] = RGCmodels.Watson.compute.coneSpacingAlongMeridianInRightEyeVisualField(-eccDegsOnHumanRetina, WatsonMeridianName, useParfor);

    % Assemble all the human temporal meridian data in one struct
    dHuman.temporalMeridianData.eccDegs = eccDegsOnHumanRetina;
    dHuman.temporalMeridianData.eccMM = eccMMOnHumanRetina;
    dHuman.temporalMeridianData.coneDensityPerSquaredMM = coneDensityMMs2;
    dHuman.temporalMeridianData.coneSpacingDegs = coneSpacingDegs;

    % Retrieve the angular eccentricity (nasal meridian)
    eccDegsOnMacaqueRetina = dMacaque.nasalMeridianData.eccDegs;
    eccMMOnMacaqueRetina = dMacaque.nasalMeridianData.eccMM;

    % Assume human and macaque retinas are equivalent at identical linear
    % ecccentricities
    eccMMOnHumanRetina = eccMMOnMacaqueRetina;

    % Derive angular eccentricity in human retina
    eccDegsOnHumanRetina = RGCmodels.Watson.convert.rhoMMsToDegs(eccMMOnHumanRetina);

    % Obtain human cone density data along the nasal meridian
    WatsonMeridianName = RGCmodels.Watson.convert.retinalMeridianNameToVisualMeridianName('nasal meridian', 'right eye');
    [coneSpacingDegs, ~, ~, coneDensityMMs2] = RGCmodels.Watson.compute.coneSpacingAlongMeridianInRightEyeVisualField(eccDegsOnHumanRetina, WatsonMeridianName, useParfor);

    % Assemble all the human nasal meridian data in one struct
    dHuman.nasalMeridianData.eccDegs = eccDegsOnHumanRetina;
    dHuman.nasalMeridianData.eccMM = eccMMOnHumanRetina;
    dHuman.nasalMeridianData.coneDensityPerSquaredMM = coneDensityMMs2;
    dHuman.nasalMeridianData.coneSpacingDegs = coneSpacingDegs;


    % Concatenate temporal and nasal meridian data
    macaqueEccMM = [flipud(dMacaque.temporalMeridianData.eccMM(:)); dMacaque.nasalMeridianData.eccMM(:)];
    macaqueEccDegs = [flipud(dMacaque.temporalMeridianData.eccDegs(:)); dMacaque.nasalMeridianData.eccDegs(:)];
    macaqueDensity  = [flipud(dMacaque.temporalMeridianData.coneDensityPerSquaredMM(:)); dMacaque.nasalMeridianData.coneDensityPerSquaredMM(:)];
    
    humanEccMM = [flipud(dHuman.temporalMeridianData.eccMM(:)); dHuman.nasalMeridianData.eccMM(:)];
    humanEccDegs = [flipud(dHuman.temporalMeridianData.eccDegs(:)); dHuman.nasalMeridianData.eccDegs(:)];
    humanDensity = [flipud(dHuman.temporalMeridianData.coneDensityPerSquaredMM(:)); dHuman.nasalMeridianData.coneDensityPerSquaredMM(:)];

    % Plot cone densities as a function of  linear eccentricity
    eccRangeMM = [-20 20];
    coneDensityRange = [1e3 3*1e5];
    eccTicksDegs = -90:30:90;
    eccTicksMM = -20:5:20;

   % eccRangeMM = [-1 1]
   % coneDensityRange = [1e4 3*1e5];
   % eccTicksDegs = -90:1:90;
   % eccTicksMM = -20:0.2:20;

    eccRangeDegs = eccRangeMM*1000/micronsPerDegreeInMacaqueRetina;

    % Prepare figure and axes
    hFig = figure(1); clf;
    ff = PublicationReadyPlotLib.figureComponents('1x1 double width figure');
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    % Plot as a function of linear eccentricity
    plot(ax,macaqueEccMM, macaqueDensity, ...
        'rs-', 'MarkerSize', ff.markerSize*0.75, 'MarkerFaceColor', [1 0.7 0.7], 'LineWidth', ff.lineWidth);
    hold(ax, 'on');
    plot(ax, humanEccMM, humanDensity, ...
        'ko-', 'MarkerSize', ff.markerSize*0.5, 'MarkerFaceColor', [0 0 0], 'LineWidth', ff.lineWidth);

    set(ax, 'YScale', 'log', 'YLim', coneDensityRange, 'YTick', [1e3 1e4 1e5 1e6]);
    set(ax, 'XLim', eccRangeMM , 'XTick', eccTicksMM);
    grid(ax, 'on'); box(ax, 'off');
    legend(ax, {'macaque (Packer)', 'human (Curcio)'});
    xlabel(ax, '\leftarrow temporal     linear ecc (mm)     nasal \rightarrow');
    ylabel(ax, 'cone density (cones/mm^2)');
    xtickangle(ax, 0)
    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);

    figureDir = fullfile(isetbioRootPath,'local',mfilename);
    if (~exist(figureDir,'dir'))
        mkdir(figureDir);
    end
    fprintf('Will save figures/videos into %s\n',figureDir);

    NicePlot.exportFigToPDF(fullfile(figureDir,'ecc_mms.pdf'),hFig,  300);

    % Plot cone densities as a function of angular eccentricity

    % Prepare figure and axes
    hFig = figure(2); clf;
    ff = PublicationReadyPlotLib.figureComponents('1x1 double width figure');
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    plot(ax,macaqueEccDegs, macaqueDensity, ...
        'rs-', 'MarkerSize', ff.markerSize*0.75, 'MarkerFaceColor', [1 0.75 0.75], 'LineWidth', ff.lineWidth);
    hold(ax, 'on');
    plot(ax, humanEccDegs, humanDensity, ...
        'ko-', 'MarkerSize', ff.markerSize*0.5, 'MarkerFaceColor', [0 0 0], 'LineWidth', ff.lineWidth);

    set(ax, 'YScale', 'log', 'YLim', coneDensityRange, 'YTick', [1e3 1e4 1e5 1e6]);
    set(ax, 'XLim', eccRangeDegs , 'XTick', eccTicksDegs);
    grid(ax, 'on'); box(ax, 'off');
    legend(ax, {'macaque (Packer)', 'human (Curcio)'});
    xlabel(ax, '\leftarrow temporal     angular ecc (degs)     nasal \rightarrow');
    ylabel(ax, 'cone density (cones/mm^2)');
    xtickangle(ax, 0)
    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);
    NicePlot.exportFigToPDF(fullfile(figureDir,'ecc_degs.pdf'),hFig,  300);

end