function [srgb,LMSphotons] = sensorDemosaicCones(sensor, method, nFrames)
% Render an image that produces a pattern of cone absorptions
%
%    [srgb,LMSphotons] = sensorDemosaicCones(sensor, [method], [nFrames])
%
% This function produces an sRGB image that creates the pattern of cone
% absorptions in the sensor object.  To create this image, 
%
%  * the cone absorption matrix is demosaicked
%  * It is then converted to an equivalent XYZ representation and 
%  * From there an sRGB image is calculated that would produce the XYZ (and
%    hence LMS) values 
%
% This function is helpful in visualizing the image one might infer given a
% particular set of cone absorptions.  It also includes the option of
% indicating that the sensor is dichromatic.  In that case the demosaicking
% from the dichromatic to sRGB image is based on the methods developed by
% HJ and by Brettel et al. and implemented in lms2lmsDichromatic.
%
%  Inputs:
%    sensor  - sensor structure with photon absorption rates computed
%    method  - algorithm to be used for interpolation. The following
%              methods are supported via griddata, applied to each
%              submosaic separately:
%                'nearest'
%                'linear' (default)
%                'natural'
%                'cubic'
%                'v4'
%              'freeman' - Freeman's method.  Median filtering applied to
%                          color plane differences after basic interpolation
%    nFrames - number of frames to be rendered, should be no larger than
%              number of frames in sensorGet(sensor, 'photons'). Default 1
%
%  Output:
%    srgb    - rendered srgb image, 4D matrix as row x col x 3 x nFrames
%    LMSphotons - The interpolated isomerizations in same format.  These
%                 are not normalized, nor promoted for dichromats. If the
%                 image is monochrome, this is returned empty.
%
%  Notes:
%    1) This function can only be used for demosaicing human cone mosaic.
%       For bayer patterns in cameras, use ISET camera modules instead
%    2) This function will detect dichromacy by checking cone densities.
%       For dichromatic observers, the rendered image will be the
%       tranformed image for trichromats. See lms2lmsDichromat for more
%       details about dichromatic color transformation
%    3) For monochrome cone mosaic, this function will just scale it a gray
%       scale image
%    4) Probably want to pass median filter size and master channel as an
%       optional parameters for method 'freeman'.  But right now they are
%       just set in that section of the code.
%
%  Examples:
%  Show the sRGB image for a trichromat
%
%    fov = 1;
%    scene = sceneCreate; scene = sceneSet(scene, 'h fov', fov);
%    oi = oiCreate('human'); oi = oiCompute(oi, scene);
%    sensor = sensorCreate('human');
%    sensor = sensorSetSizeToFOV(sensor, fov, scene, oi);
%    sensor = sensorCompute(sensor, oi);
%    srgb = sensorDemosaicCones(sensor, 'linear');
%    vcNewGraphWin; imshow(srgb);
%
%  Shoe the sRGB image for a protanope, with no noise
%    cP = coneCreate;
%    sensor = sensorCreate('human',cP);
%    sensor = sensorSetSizeToFOV(sensor, fov, scene, oi);
%    sensor = sensorCompute(sensor, oi);
%    srgb = sensorDemosaicCones(sensor, 'linear');
%    vcNewGraphWin; imshow(srgb);

%    sensor = sensorComputeNoiseFree(sensor, oi);
%    srgb   = sensorDemosaicCones(sensor, 'linear');
%    vcNewGraphWin; imshow(srgb);
%
%  See also:
%    lms2lmsDichromat, sensorGet, sensorPlot
%
% (HJ) ISETBIO TEAM, 2015
%
% 4/16/15  dhb  Added Freeman method.

%% Check inputs
if notDefined('sensor'), error('sensor required'); end
if notDefined('method'), method = 'linear'; end
if notDefined('nFrames'), nFrames = 1; end
if ~sensorCheckHuman(sensor), error('only human sensor supported'); end

%% Set up for various methods
switch (method)
    case { 'nearest' 'linear' 'natural' 'cubic' 'v4'}
        gdmethod = method;
    case {'freeman'}
        gdmethod = 'linear';
    otherwise
        error('Unknown method passed');
end

p = sensorGet(sensor, 'photons');
if isempty(p), error('photon absorption not computed'); end
if nFrames > size(p, 3), error('nFrames exceeds photons frames'); end

p  = p(:, :, 1:nFrames); % use first nFrames only
sz = sensorGet(sensor, 'size'); % cone mosaic size

%% Check mosaic type (trichromats, dichromats, monochrome)
%  get cone spatial densities
density  = sensorGet(sensor, 'human cone densities'); % KLMS proportions

%  check number of non-zero entries
cbType = find(~density(2:4));
if length(cbType) > 2, error('Bad cone density: %.2f', density); end

%% Interpolate
if length(cbType) == 2 % monochrome case
    % scale photon absorptions to get the rendered image
    srgb = p / max(p(:));
    srgb = reshape(srgb, [sz 1 nFrames]);
    srgb = repmat(srgb, [1 1 3 1]);
    LMSphotons = [];
else % trichromatic or dichromatic case
    coneType = sensorGet(sensor, 'cone type');
    
    % allocate space
    srgb = zeros([sz 3 nFrames]);
    LMSphotons  = zeros([sz 3 nFrames]);
    
    % Find scaling so that conversion to sRGB is correct. That conversion
    % assumes that LMS spectral sensitivities are scaled to peak of 1 in
    % energy units. Here we have real quantal efficiences, so that our LMS
    % values are scaled differently relative to one another. This will
    % produce problems if we just use the standard transformation matrix.
    %
    % Getting this totally right is a bit involved, particularly if we want
    % to get the XYZ values in standard units of Y being in cd/m2 or some
    % such. Here we'll only fuss a little, and get the relative scaling
    % right.
    %
    % Note that we need to convert to energy units to make the relative
    % scaling right.
    sensorQE = sensorGet(sensor, 'spectral qe');
    wave = sensorGet(sensor,'wave');
    sensorEnergyUnits = EnergyToQuanta(wave(:),sensorQE);
    scale = max(sensorEnergyUnits); 
    scale = scale/max(scale);
    scale = reshape(scale(2:4), [1 1 3]);
    
    % Set up interpolation mesh
    [xq, yq] = meshgrid(1:sz(2), 1:sz(1));

    % loop over frames
    for curFrame = 1 : nFrames 
        curP = p(:,:,curFrame);
        % interpolate for L, M and S, which are indexed as 2:4 in isetbio.
        for ii = 2 : 4 
            indx = find(coneType == ii); 
            if isempty(indx), LMSphotons(:,:,ii-1) = 0; continue; end
            
            % ind2sub works in i, j so that the return arguments should be
            % interpreted as [y, x] and not the other way around
            [y, x] = ind2sub(sz, indx);
            
            % Let griddata do the work.
            LMSphotons(:,:,ii-1,curFrame) = griddata(x,y,curP(indx), xq, yq, gdmethod);
        end
        
        % Adjust interpolation if it's freeman's method
        if (strcmp(method,'freeman'))
            masterSensor = 1;
            filterSize = 5;
            freemanInterpLMS = LMSphotons(:,:,:,curFrame);
            tempMaster = LMSphotons(:,:,masterSensor,curFrame);
            
            % Leave the masterSensor plane along, do median filtering on 
            % difference images for the other two.  For each plane, use
            % the median filtered version at locations where there isn't a
            % sensor of that plane.
            %
            % This may not be the exact algorithm in the Freeman patent.
            % My memory is that Bill suggested this version to me as one
            % that was easy to implement and worked just as well.
            %
            % Although the Freeman algorithm works pretty well for the
            % Bayer mosaic, it only appears to help a little for human cone
            % layouts.  It is possible this implementation is buggy.
            for n = 1:3
                if (n ~= masterSensor)
                    diffImage = LMSphotons(:,:,n,curFrame) - tempMaster;
                    diffImage = medfilt2(diffImage,[filterSize filterSize]);
                    index = find(coneType ~= n+1);
                    temp = LMSphotons(:,:,n,curFrame);
                    temp(index) = diffImage(index) + tempMaster(index);
                    freemanInterpLMS(:,:,n) = temp;
                end
            end
            
            % Pop this back in
            LMSPhotons(:,:,:,curFrame) = freemanInterpLMS;
        end
        
        % scale LMS to match stockman normalized spectral qe
        % we need to do this because all xyz2lms and colorblind transform
        % routine are based on stockman normalized spectral qe
        LMS = bsxfun(@rdivide, LMSphotons(:,:,:,curFrame), scale);
        
        % colorblind transformation
        LMS = lms2lmsDichromat(LMS, cbType, 'linear');
        
        % transform back to srgb
        srgb(:,:,:,curFrame) = lms2srgb(LMS);
    end
end

end