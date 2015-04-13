function srgb = sensorDemosaicCones(sensor, method, nFrames)
%% Demosaic and render human cone photon absorptions
%    srgb = sensorDemosaicCones(sensor, [method], [nFrames])
%
%  Inputs:
%    sensor  - sensor structure with photon absorption rates computed
%    method  - algorithm to be used for interpolation. Now support
%              'nearest', 'linear' (default), 'natural', 'cubic', 'v4'
%    nFrames - number of frames to be rendered, should be no larger than
%              number of frames in sensorGet(sensor, 'photons'). Default 1
%
%  Output:
%    srgb    - rendered srgb image, 4D matrix as row x col x 3 x nFrames
%
%  Notes:
%    1) this function can only be used for demosaicing human cone mosaic.
%       For bayer patterns in cameras, use ISET camera modules instead
%    2) this function will detect dichromacy by checking cone densities.
%       For dichromatic observers, the rendered image will be the
%       tranformed image for trichromats. See lms2lmsDichromat for more
%       details about dichromatic color transformation
%    3) For monochrome cone mosaic, this function will just scale it a gray
%       scale image
%
%  Examples:
%    fov = 1;
%    scene = sceneCreate; scene = sceneSet(scene, 'h fov', fov);
%    oi = oiCreate('human'); oi = oiCompute(oi, scene);
%    sensor = sensorCreate('human');
%    sensor = sensorSetSizeToFOV(sensor, fov, scene, oi);
%    sensor = sensorCompute(sensor, oi);
%    srgb = sensorDemosaicCones(sensor, 'linear');
%    vcNewGraphWin; imshow(srgb);
%
%    sensor = sensorComputeNoiseFree(sensor, oi);
%    srgb   = sensorDemosaicCones(sensor, 'linear');
%    vcNewGraphWin; imshow(srgb);
%
%  See also:
%    lms2lmsDichromat, sensorGet, plotSensor
%
% (HJ) ISETBIO TEAM, 2015

%% Check inputs
if notDefined('sensor'), error('sensor required'); end
if notDefined('method'), method = 'linear'; end
if notDefined('nFrames'), nFrames = 1; end
if ~sensorCheckHuman(sensor), error('only human sensor supported'); end

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
else % trichromatic or dichromatic case
    coneType = sensorGet(sensor, 'cone type');
    
    % allocate spaces
    srgb = zeros([sz 3 nFrames]);
    LMS  = zeros([sz 3]);
    
    % define const parameters
    [yq, xq] = meshgrid(1:sz(2), 1:sz(1));
    scale = max(sensorGet(sensor, 'spectral qe'));
    scale = reshape(scale(2:4), [1 1 3]);
    
    for curFrame = 1 : nFrames % loop over frames
        curP = p(:,:,curFrame);
        for ii = 2 : 4 % interpolate for L, M and S
            indx = find(coneType == ii); 
            if isempty(indx), LMS(:,:,ii-1) = 0; continue; end
            [x, y] = ind2sub(sz, indx);
            LMS(:,:,ii-1) = griddata(x,y,curP(indx), xq, yq, method);
        end
        % scale LMS to match stockman normalized spectral qe
        % we need to do this because all xyz2lms and colorblind transform
        % routine are based on stockman normalized spectral qe
        LMS = bsxfun(@rdivide, LMS, scale);
        
        % colorblind transformation
        LMS = lms2lmsDichromat(LMS, cbType, 'linear');
        
        % transform back to srgb
        srgb(:,:,:,curFrame) = lms2srgb(LMS);
    end
end

end