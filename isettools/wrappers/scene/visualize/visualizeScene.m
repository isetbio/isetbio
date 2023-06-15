function visualizeScene(scene, varargin)
    p = inputParser;
    p.addParameter('presentationDisplay', [], @(x)(isempty(x)||isstruct(x)));
    p.addParameter('displayContrastProfiles', false, @islogical);
    p.addParameter('displayRadianceMaps', true, @islogical);
    p.addParameter('spatialSupportInDegs', false, @islogical);
    p.addParameter('avoidAutomaticRGBscaling', false, @islogical);
    p.addParameter('roiRectDegs', struct(), @isstruct);
    p.addParameter('crossHairsAtOrigin', false, @islogical);
    p.addParameter('axesHandle', []);
    p.addParameter('noTitle', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.addParameter('noYTicks', false, @islogical);
    % Parse input
    p.parse(varargin{:});
    
    presentationDisplay = p.Results.presentationDisplay;
    displayContrastProfiles = p.Results.displayContrastProfiles;
    displayRadianceMaps = p.Results.displayRadianceMaps;
    spatialSupportInDegs = p.Results.spatialSupportInDegs;
    roiRectDegs = p.Results.roiRectDegs;
    crossHairsAtOrigin = p.Results.crossHairsAtOrigin;
    avoidAutomaticRGBscaling = p.Results.avoidAutomaticRGBscaling;
    

    % retrieve the spatial support of the scene(in millimeters)
    spatialSupportMilliMeters = sceneGet(scene, 'spatial support', 'mm');

    % retrieve the XYZ tristimulus components of the scene
    XYZmap = sceneGet(scene, 'xyz');
    % compute the xy-chroma and the luminance maps
    luminanceMap = squeeze(XYZmap(:,:,2));
    xMap = squeeze(XYZmap(:,:,1))./sum(XYZmap,3);
    yMap = squeeze(XYZmap(:,:,2))./sum(XYZmap,3);
    % Compute mean luminance and mean chromaticity
    meanLuminance = mean(luminanceMap(:));
    meanChromaticity = [mean(xMap(:)) mean(yMap(:))];

    % retrieve the sRGB components of the scene (just for visualization)
    if (isempty(presentationDisplay))
        sceneRGBsettings = sceneGet(scene, 'rgb image');
    else
        % If the presentationDisplay that was used to generate this scene
        % is passed we can do a much better job
        % computing the sceneRGBsettings from the sceneXYZ
        displayLinearRGBToXYZ = displayGet(presentationDisplay, 'rgb2xyz');
        displayXYZToLinearRGB = inv(displayLinearRGBToXYZ);
        sceneRGBsettings = imageLinearTransform(XYZmap, displayXYZToLinearRGB);
    end

    if (spatialSupportInDegs)
        viewingDistance = sceneGet(scene, 'distance');
        spatialSupportDegs = 2 * atand(spatialSupportMilliMeters/1e3/2/viewingDistance);
    
        spatialSupport = spatialSupportDegs;
        spatialSupportUnits = 'degs';
    else
        spatialSupport = spatialSupportMilliMeters;
        spatialSupportUnits = 'mm';
    end
    
    % visualize the scene as RGB
    ax = visualizeSceneRGB(spatialSupport, spatialSupportUnits, sceneRGBsettings, ...
        meanLuminance, meanChromaticity, sceneGet(scene, 'name'), ...
        'avoidAutomaticRGBscaling', avoidAutomaticRGBscaling, ...
        'axesHandle', p.Results.axesHandle, ...
        'noTitle', p.Results.noTitle, ...
        'noYLabel', p.Results.noYLabel, ...
        'noYTicks', p.Results.noYTicks);
    
    
    % Add ROI
    if (~isempty(fieldnames(roiRectDegs))) && (spatialSupportInDegs)
        xRect = [...
            roiRectDegs.xo + roiRectDegs.width/2*(-1), ...
            roiRectDegs.xo + roiRectDegs.width/2*(-1), ...
            roiRectDegs.xo + roiRectDegs.width/2, ...
            roiRectDegs.xo + roiRectDegs.width/2, ...
            roiRectDegs.xo + roiRectDegs.width/2*(-1) ...
            ];
        yRect = [...
            roiRectDegs.yo + roiRectDegs.height/2*(-1), ...
            roiRectDegs.yo + roiRectDegs.height/2, ...
            roiRectDegs.yo + roiRectDegs.height/2, ...
            roiRectDegs.yo + roiRectDegs.height/2*(-1), ...
            roiRectDegs.yo + roiRectDegs.height/2*(-1) ...
            ];
        hold(ax, 'on');
        plot(ax, xRect,yRect, 'r-', 'LineWidth', 1.5);
    end
    
    if (crossHairsAtOrigin)
        spatialSupportX = squeeze(spatialSupport(1,:,1));
        spatialSupportY = squeeze(spatialSupport(:,1,2));
        hold(ax, 'on');
        plot(ax,[spatialSupportX(1) spatialSupportX(end)], [0 0], 'k-');
        plot(ax,[0 0], [spatialSupportY(1) spatialSupportY(end)],  'k-');
    end
    
    if (displayRadianceMaps)
        % retrieve the radiance of the scene as emitted photon rate 
        % (photons/pixel/sec/nm)
        scenePhotonRate = sceneGet(scene, 'photons');
        % Retrieve wavelength support
        wavelengthSupport = sceneGet(scene, 'wave');
        wavelengthBandsToVisualize = 400:30:700;
    
        % Visualize radiance maps
        visualizeSceneRadiance(spatialSupportMilliMeters, 'mm', ...
            scenePhotonRate, wavelengthSupport, wavelengthBandsToVisualize);
    end
    
    if (displayContrastProfiles)
        visualizeSceneRadiance(spatialSupportMilliMeters, 'mm', ...
        scenePhotonRate, wavelengthSupport, wavelengthBandsToVisualize, ...
        'contrastProfilesOnly', true); 
    end

    drawnow;
end

