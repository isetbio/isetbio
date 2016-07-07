% t_osHyperSpectralSceneEyeScan
%
% A tutorial script illustrating usage of the outer segment object in
% computing cone responses to eye scanning of hyperspectral images. 
%
% Dependencies: RemoteDataToolbox for fetching hyperspectral scenes
%
% 1/2016        NPC   Created.
% 2/24/2016     NPC   Added option to generate video (good for presentations)

function t_osHyperSpectralSceneEyeScan
 
    % reset
    ieInit; close all;

    % Set up remote data toolbox client
    client = RdtClient('isetbio');

    % Spacify images - Updataed/HJ/BW
    imageSources = {...
            {'manchester_database/2002', 'scene1_2002_iset'} ...
    };

    % simulation time step. same for eye movements and for sensor, outersegment
    timeStepInMilliseconds = 0.02;

    sensorParams = struct(...
        'coneApertureInMicrons', 3.0, ...        % custom cone aperture
        'LMSdensities', [0.6 0.4 0.1], ...       % custom percentages of L,M and S cones
        'spatialGrid', [20 20], ...              % generate mosaic of 20 x 20 cones
        'samplingIntervalInMilliseconds', timeStepInMilliseconds, ...
        'integrationTimeInMilliseconds', 50, ...
        'randomSeed', 1552784, ...
        'eyeMovementScanningParams', struct(...
            'samplingIntervalInMilliseconds', timeStepInMilliseconds, ...
            'fixationDurationInMilliseconds', 300, ...
            'numberOfFixations', 3 ...
        ) ...
    );

    % go through all the images
    for imageIndex = 1:numel(imageSources)
        % retrieve scene name
        imsource = imageSources{imageIndex};
        hProgress = waitbar(0.1, sprintf('Fetching scene named ''%s'' (''%s'') from isetbio''s remote repository. Please wait ...', imsource{2}, imsource{1}));

        % download image - Updated HJ/BW
        client.crp(sprintf('/resources/scenes/hyperspectral/%s', imsource{1}));
        data = client.readArtifact(imsource{2}, 'type', 'mat');
        scene = sceneFromBasis(data);
        clear data;
        % scene = artifactData.scene; clear 'artifactData'; clear 'artifactInfo';

        % adjust mean luminance to avoid outer segment photocurrent saturation
        forcedMeanLuminanceInCdPerM2 = min([sceneGet(scene, 'mean luminance') 600]);
        scene = sceneAdjustLuminance(scene, forcedMeanLuminanceInCdPerM2);
        
        % Show scene
        waitbar(0.25, hProgress,'Displaying scene ...'); figure(hProgress);
        %vcAddAndSelectObject(scene); sceneWindow;

        % Compute optical image with human optics
        waitbar(0.3, hProgress,'Computing optical image ...'); figure(hProgress);
        oi = oiCreate('human');
        oi = oiCompute(oi, scene);

        % Show optical image
        %vcAddAndSelectObject(oi); oiWindow;

        % create custom human sensor
        waitbar(0.5, hProgress,'Computing photoreceptor isomerizations ...'); figure(hProgress);
        sensor = sensorCreate('human');
        sensor = customizeSensor(sensor, sensorParams, oi);

        % compute rate of isomerized photons
        sensor = coneAbsorptions(sensor, oi);

        % compute outer segment response using the linear adaptation model
        waitbar(0.7, hProgress,'Computing linear outer segment response ...'); figure(hProgress);
        osL = osLinear(); 
        osL = osSet(osL, 'noiseFlag', 1);
        osL = osCompute(osL, sensor);
        % display osWindow for interactive viewing of outer segment responses
        osLwindow = osWindow(1001+imageIndex, 'linear outer segment', 'horizontalLayout', osL, sensor, oi, scene);

        % compute outer segment response using the biophysically-based adaptation model
        waitbar(0.9, hProgress,'Computing biophysically-based outer segment response ...'); figure(hProgress);
        osB = osBioPhys();
        osB = osSet(osB, 'noiseFlag', 1);
        osB = osCompute(osB, sensor);
        
        % display osWindow for interactive viewing of outer segment responses
        osBwindow = osWindow(1002+imageIndex, 'biophys-based outer segment', 'horizontalLayout', osB, sensor, oi, scene);
        close(hProgress);
         
        % video generation: good for presentations
        genVideo = input('Generate video of the linear outer segment simulation ? [y/n] : ', 's');
        if (strcmp(genVideo, 'y'))
            % Generate video stream, one frame/msec of simulation
            videoFrameStepInMilliseconds = 1.0;
            videoFilename = sprintf('LinearOSSimulation_%s.m4v', imsource{2});
            osLwindow.generateVideo(videoFrameStepInMilliseconds, videoFilename);
        end
        
        genVideo = input('Generate video of the biophys-based outer segment simulation ? [y/n] : ', 's');
        if (strcmp(genVideo, 'y'))
            % Generate video stream, one frame/msec of simulation
            videoFrameStepInMilliseconds = 1.0;
            videoFilename = sprintf('BiophysOSSimulation_%s.m4v', imsource{2});
            osBwindow.generateVideo(videoFrameStepInMilliseconds, videoFilename);
        end
        
        % test resizing
        testFigureResize = false;
        if (testFigureResize)
            while(1)
                desiredWindowSize = input('Enter new figure size in pixels: ');
                osLwindow.resizeWidow(desiredWindowSize);
            end
        end
        
    end
end


function sensor = customizeSensor(sensor, sensorParams, opticalImage)
    
    if (isempty(sensorParams.randomSeed))
       rng('shuffle');   % produce different random numbers
    else
       rng(sensorParams.randomSeed);
    end
    
    eyeMovementScanningParams = sensorParams.eyeMovementScanningParams;
    
    % custom aperture
    pixel  = sensorGet(sensor,'pixel');
    pixel  = pixelSet(pixel, 'size', [1.0 1.0]*sensorParams.coneApertureInMicrons*1e-6);  % specified in meters);
    sensor = sensorSet(sensor, 'pixel', pixel);
    
    % custom LMS densities
    coneMosaic = coneCreate();
    coneMosaic = coneSet(coneMosaic, ...
        'spatial density', [0.0 ...
                           sensorParams.LMSdensities(1) ...
                           sensorParams.LMSdensities(2) ...
                           sensorParams.LMSdensities(3)] );
    sensor = sensorCreateConeMosaic(sensor,coneMosaic);
        
    % sensor wavelength sampling must match that of opticalimage
    sensor = sensorSet(sensor, 'wavelength', oiGet(opticalImage, 'wavelength'));
     
    % no noise on sensor
    sensor = sensorSet(sensor,'noise flag', 0);
    
    % custom size
    sensor = sensorSet(sensor, 'size', sensorParams.spatialGrid);

    % custom time interval
    sensor = sensorSet(sensor, 'time interval', sensorParams.samplingIntervalInMilliseconds/1000.0);
    
    % custom integration time
    sensor = sensorSet(sensor, 'integration time', sensorParams.integrationTimeInMilliseconds/1000.0);
    
    % custom eye movement
    eyeMovement = emCreate();
    
    % custom sample time
    eyeMovement  = emSet(eyeMovement, 'sample time', eyeMovementScanningParams.samplingIntervalInMilliseconds/1000.0);        
    
    % attach eyeMovement to the sensor
    sensor = sensorSet(sensor,'eyemove', eyeMovement);
            
    % generate the fixation eye movement sequence
    eyeMovementsNum = eyeMovementScanningParams.numberOfFixations * round(eyeMovementScanningParams.fixationDurationInMilliseconds / eyeMovementScanningParams.samplingIntervalInMilliseconds);
    eyeMovementPositions = zeros(eyeMovementsNum,2);
    sensor = sensorSet(sensor,'positions', eyeMovementPositions);
    sensor = emGenSequence(sensor);
    
    % add saccadic targets
    saccadicTargetPos = generateSaccadicTargets(eyeMovementScanningParams, 'random'); 
    eyeMovementPositions = sensorGet(sensor,'positions', eyeMovementPositions);
    for eyeMovementIndex = 1:eyeMovementsNum
        kk = 1+floor((eyeMovementIndex-1)/round(eyeMovementScanningParams.fixationDurationInMilliseconds / eyeMovementScanningParams.samplingIntervalInMilliseconds));
        eyeMovementPositions(eyeMovementIndex,:) = eyeMovementPositions(eyeMovementIndex,:) + saccadicTargetPos(kk,:);
    end
    sensor = sensorSet(sensor,'positions', eyeMovementPositions);
end

function saccadicTargetPos = generateSaccadicTargets(eyeMovementScanningParams, mode) 

    saccadicTargetPos = zeros(eyeMovementScanningParams.numberOfFixations,2);
    
    % random targets
    if (strcmp(mode, 'random'))
        saccadicTargetPos = round(randn(eyeMovementScanningParams.numberOfFixations,2)*100);
    else
        % oscillate between three positions of interest - useful for debuging
        for k = 1:eyeMovementScanningParams.numberOfFixations
            if (mod(k-1,6) < 2)
                saccadicTargetPos(k,:) = [-850 -390]/3;  
            elseif (mod(k-1,6) < 4)
                saccadicTargetPos(k,:) = [-170 515]/3;
            else
                saccadicTargetPos(k,:) = [-105 505]/3; 
            end
        end
    end
end

    
        
    