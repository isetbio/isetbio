% Demo usage of rodIntrusionAdjustedConeAperture flag in new @cMosaic object
%
% Description:
%    Shows how to generate @cMosaics in which the cone aperture takes into
%    account rod intrusion as a function of eccentricity. Also compares the
%    generated mosaics with those of Curcio 1990. (Fig 3).
%

% History:
%    07/15/22  NPC  ISETBIO Team, Copyright 2021 Wrote it.

function t_cMosaicRodIntrusion()

    % Load Figure 3 of Curcio 1990 
    resourcesDir = strrep(isetRootPath, 'isettools', 'tutorials/cones/resources');
    load(fullfile(resourcesDir, 'CurcioConesRods.mat'), 'CurcioConesRods');

    hFig = figure(1);clf;
    set(hFig, 'Position', [70 70 1280 890]);

    % Plot it
    ax = subplot('Position', [0.01 0.01 0.48 0.98]);
    imagesc(ax,CurcioConesRods);
    axis(ax,'image');
    colormap(gray);
    set(ax, 'XTick', [], 'YTick', []);

    % Curcio 1990, Figure 3, mosaic eccentricities
    eccMM = [1.35 5.0 8.0];

    % Transform to degrees
    eccDegs = RGCmodels.Watson.convert.rhoMMsToDegs(eccMM);   

    % Positions for the various subplots
    sv = NicePlot.getSubPlotPosVectors(...
           'colsNum', 2, ...
           'rowsNum', 4, ...
           'heightMargin',  0.02, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.04, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.06, ...
           'topMargin',      0.09); 

    whichEye = 'left eye';
    for positionIndex = 1:numel(eccDegs)

        % Generate nasal mosaic at specified eccentricity
        [horizontalEccDegs, verticalEccDegs] = cMosaic.eccentricitiesForRetinaMeridianInEye(...
            eccDegs(positionIndex), 'nasal meridian', whichEye);

        cMosaicNasalMeridian{positionIndex} = cMosaic(...
            'eccentricityDegs', [horizontalEccDegs verticalEccDegs], ...
            'sizeDegs', [0.5 0.5], ...
            'coneDensities', [1.0 0 0 0], ...
            'whichEye', whichEye, ...
            'rodIntrusionAdjustedConeAperture', true, ...
            'coneApertureModifiers', struct('smoothLocalVariations', true));

        % Generate temporal mosaic at specified eccentricity
        [horizontalEccDegs, verticalEccDegs] = cMosaic.eccentricitiesForRetinaMeridianInEye(...
            eccDegs(positionIndex), 'temporal meridian', whichEye);

        cMosaicTemporalMeridian{positionIndex} = cMosaic(...
            'eccentricityDegs', [horizontalEccDegs verticalEccDegs], ...
            'sizeDegs', [0.5 0.5], ...
            'coneDensities', [1.0 0 0 0], ... 
            'whichEye', whichEye, ...
            'rodIntrusionAdjustedConeAperture', true, ...
            'coneApertureModifiers', struct('smoothLocalVariations', true));

        % XY limits for the nasal mosaic
        xyLimitsNasal(positionIndex,:) = [...
            cMosaicNasalMeridian{positionIndex}.eccentricityMicrons(1)-19 ...
            cMosaicNasalMeridian{positionIndex}.eccentricityMicrons(1)+19 ...
            cMosaicNasalMeridian{positionIndex}.eccentricityMicrons(2)-13 ...
            cMosaicNasalMeridian{positionIndex}.eccentricityMicrons(2)+13]';


        % Visualize a patch of the nasal mosaic
        axisPos = sv(positionIndex,1).v;
        axisPos(3) = 0.51*axisPos(3);
        axisPos(1) = 0.48 + axisPos(1);
        ax = subplot('Position', axisPos);
        
        plotTitle = [];
        if (positionIndex == 1)
            plotTitle = 'nasal retina';
        end

        cMosaicNasalMeridian{positionIndex}.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax, ...
            'domain', 'microns', ...
            'noXLabel', (positionIndex < size(eccDegs,1)), ...
            'domainVisualizationLimits', xyLimitsNasal(positionIndex,:), ...
            'domainVisualizationTicks', ...
                    struct(...
                        'x', round(cMosaicNasalMeridian{positionIndex}.eccentricityMicrons(1)+(-20:10:20)), ...
                        'y', round(cMosaicNasalMeridian{positionIndex}.eccentricityMicrons(2)+(-20:10:20))), ...
            'fontSize', 12, ...
            'plotTitle', plotTitle);

        % XY limits for the temporal mosaic
        xyLimitsTemporal(positionIndex,:) = [...
            cMosaicTemporalMeridian{positionIndex}.eccentricityMicrons(1)-19 ...
            cMosaicTemporalMeridian{positionIndex}.eccentricityMicrons(1)+19 ...
            cMosaicTemporalMeridian{positionIndex}.eccentricityMicrons(2)-13 ...
            cMosaicTemporalMeridian{positionIndex}.eccentricityMicrons(2)+13]';

        % Visualize a patch of the temporal mosaic
        plotTitle = [];
        if (positionIndex == 1)
            plotTitle = 'temporal retina';
        end
        xyLimitsTemporal(positionIndex,:)
        axisPos = sv(positionIndex,2).v;
        axisPos(3) = 0.51*axisPos(3);
        axisPos(1) = 0.75+(axisPos(1)-0.5);
        ax = subplot('Position', axisPos);
        cMosaicTemporalMeridian{positionIndex}.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax, ...
            'domain', 'microns', ...
            'domainVisualizationLimits', xyLimitsTemporal(positionIndex,:), ...
            'domainVisualizationTicks', ...
                    struct(...
                        'x', round(cMosaicTemporalMeridian{positionIndex}.eccentricityMicrons(1)+(-20:10:20)), ...
                        'y', round(cMosaicTemporalMeridian{positionIndex}.eccentricityMicrons(2)+(-20:10:20))), ...
            'fontSize', 12, ...
            'noXLabel', (positionIndex < size(eccDegs,1)), ...
            'noYLabel', true, ...
            'plotTitle', plotTitle);

    end

    % Show mosaic activations to a uniform field
    hFig = figure(3);clf;
    set(hFig, 'Position', [70 70 720 950]);

    uniformFieldOI = generateTestOI();

    for positionIndex = 1:numel(eccDegs)

        % Compute response to test stimulus
        activation = cMosaicNasalMeridian{positionIndex}.compute(uniformFieldOI, 'opticalImagePositionDegs', 'mosaic-centered');
    

        ax = subplot('Position', sv(positionIndex,1).v);
        plotTitle = [];
        if (positionIndex == 1)
            plotTitle = 'nasal meridian';
        end

        cMosaicNasalMeridian{positionIndex}.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax, ...
            'domain', 'microns', ...
            'domainVisualizationLimits', xyLimitsNasal(positionIndex,:), ...
            'domainVisualizationTicks', ...
                    struct(...
                        'x', round(cMosaicNasalMeridian{positionIndex}.eccentricityMicrons(1)+(-20:10:20)), ...
                        'y', round(cMosaicNasalMeridian{positionIndex}.eccentricityMicrons(2)+(-20:10:20))), ...
            'fontSize', 12, ...
            'activation', activation, ...
            'activationRange', [40 400], ...
            'noXLabel', (positionIndex < size(eccDegs,1)), ...
            'plotTitle', plotTitle);


        % Compute response to test stimulus
        activation = cMosaicTemporalMeridian{positionIndex}.compute(uniformFieldOI, 'opticalImagePositionDegs', 'mosaic-centered');

   
        plotTitle = [];
        if (positionIndex == 1)
            plotTitle = 'temporal meridian';
        end
        ax = subplot('Position', sv(positionIndex,2).v);
    
        cMosaicTemporalMeridian{positionIndex}.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax, ...
            'domain', 'microns', ...
            'activation', activation, ...
            'activationRange', [40 400], ...
            'domainVisualizationLimits', xyLimitsTemporal(positionIndex,:), ...
            'domainVisualizationTicks', ...
                    struct(...
                        'x', round(cMosaicTemporalMeridian{positionIndex}.eccentricityMicrons(1)+(-20:10:20)), ...
                        'y', round(cMosaicTemporalMeridian{positionIndex}.eccentricityMicrons(2)+(-20:10:20))), ...
            'fontSize', 12, ...
            'noXLabel', (positionIndex < size(eccDegs,1)), ...
            'noYLabel', true, ...
            'plotTitle', plotTitle);

end

% Generate a mosaic with a custom factor for cone aperture adjustment due to rod intrusion 
cMosaicTest = cMosaic(...
        'eccentricityDegs', [10 0], ...
        'sizeDegs', [0.5 0.5], ...
        'coneDensities', [1.0 0 0 0], ...
        'whichEye', whichEye, ...
        'rodIntrusionAdjustedConeAperture', 0.6, ...
        'coneApertureModifiers', struct('smoothLocalVariations', true));
cMosaicTest.visualize()

end

function oi = generateTestOI()
    pixelsNum = 512;
    fovDegs = 2.0;
    stimFreqCyclesPerDeg = 0;
    
    parms.freq = stimFreqCyclesPerDeg*fovDegs;
    parms.contrast = 1;
    parms.ph = 0;
    parms.ang = 0;
    parms.row = pixelsNum;
    parms.col = pixelsNum;
    parms.GaborFlag = 0;
    scene = sceneCreate('harmonic', parms);
    scene = sceneSet(scene, 'fov', fovDegs);
    
    %% Compute the optical image
    oi = oiCreate;
    oi = oiCompute(scene, oi);
end