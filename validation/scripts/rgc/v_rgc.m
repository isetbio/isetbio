function varargout = v_rgc(varargin)
%
% Validate the RGC object.


    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

%% Init
ieInit;



%% Compare spatial response to spatial receptive field and
%%% temporal impulse response to temporal response with impulse stimulus

%%% build scene
params.image_size = 63;
params.meanLuminance = 1;
params.nsteps = 10;
params.fov = 0.8;

% extent = round(size(rgc1.mosaic{1}.sRFcenter{1,1},1)/rgc1.mosaic{1}.rfDiameter)
% centerCoord = 21;%round(extent*rgc1.mosaic{1}.rfDiameter/2)
% [scene, display] = sceneHorwitzHassWhiteNoise(params);

% sceneRGB = sceneHorwitzHassWhiteNoiseRGB(params);
% sceneRGB = sceneHorwitzHassBarRGB(params);
sceneRGB = zeros(params.image_size, params.image_size, params.nsteps, 3);
sceneRGB(16,16,1,:) = [1 1 1];
% sceneRGB(38,38,1,:) = [1 1 1];
scene = sceneFromFile(zeros(params.image_size,params.image_size), 'rgb', params.meanLuminance);

oi  = oiCreate('wvf human');
sensor = sensorCreate('human');

sensor = sensorSet(sensor, 'integration time', .01);

identityOS = osCreate('identity');
identityOS = osSet(identityOS, 'rgbData', sceneRGB);

%%% build rgc
% rgc1 = rgcLinear(scene, sensor, osIdentity, 'right', 3.75, 180);
% rgc1 = rgcLNP(scene, sensor, osIdentity, 'right', 3.75, 180);
rgc1 = rgcGLM(scene, sensor, identityOS, 'right', 3.0, 180);
% rgc1 = rgcSubunit(scene, sensor, identityOS, 'right', 3.0, 180);

%%% compute rgc
rgc1 = rgcCompute(rgc1, identityOS);

% In order for validation to work, impulse stimulus must not be set to zero
% mean, RFs must have zero noise in center position.
rf1 = (rgc1.mosaic{1}.sRFcenter{1,1} - rgc1.mosaic{1}.sRFsurround{1,1}); rf1 = rf1./max(abs(rf1(:)));
rf2 = (rgc1.mosaic{1}.linearResponse{1,1,2}(:,:,1,1)); rf2 = rf2./max(abs(rf2(:)));
% figure; subplot(121); imagesc(rf1); subplot(122); imagesc(rf2)
diff1 = rf1 + rf2;  
max(abs(diff1(:)));

% In order for temporal validation, need to use the full length of the
% temporal response, not just the center/same. This needs to be adjusted in
% fullConvolve.m.
ir1 = (rgc1.mosaic{1}.tCenter{1,1}); ir1 = ir1./max(abs(ir1(:)));
ir2 = (rgc1.mosaic{1}.linearResponse{1,1,1}(:,:,1)); ir2 = ir2./max(abs(ir2(:)));
diff2 = ir1(:)' - ir2(1:length(ir1));  
max(abs(diff2(:)));

%%
tolerance = 1e-14;
UnitTest.assertIsZero(max(abs(diff1(:))),'Comparison for spatial RF and spatial impulse response',tolerance);
UnitTest.validationData('spatial RF',rf1);
UnitTest.validationData('spatial impulse response',rf2);

tolerance = 1e-14;
UnitTest.assertIsZero(max(abs(diff2(:))),'Comparison for temporal IR and temporal response to impulse',tolerance);
UnitTest.validationData('temporal IR',rf1);
UnitTest.validationData('temporal response to impulse',rf2);

end

